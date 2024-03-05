#include <iostream>
#include <Eigen/Dense>
#include "ezOptionParser.hpp"
#include <omp.h>
#include <sys/stat.h>

#include "utils.h"

void pod(ez::ezOptionParser &opt)
{
  std::cout << "Starting POD routine " << std::endl ;

  Parameters params(opt) ;

  omp_set_num_threads(params.m_threadsSize) ;

  // GENERATING TIME STRING
  std::vector<std::string> t(read_timefile(params.m_timesFileName)) ;

  long timesSize(t.size()) ; // Determine number of snapshots from time entry list

  /* Check whether the requested number of modes to write is larger than
  the number of snapshots. If so, set the number of modes to the number of
  snapshots. */
  if (params.m_podSize > timesSize)
  {
    std::cout << "Modes to write exceed available snapshots. Adjusted. " << std::endl ;
    params.m_podSize = timesSize ;
  }

  // READING INPUT FILES

  /* Build a string array with all file names to be read into the matrix */
  std::vector<std::string> pcfs ;
  for (std::vector<std::string>::iterator it = t.begin(); it != t.end(); ++it)
  {
    std::string name_temp = params.m_inputDirName + "/" + *it + "/" + params.m_dataFileName ;
    pcfs.push_back(name_temp);
  }

  MatrixXd m ;
  double start(omp_get_wtime()) ;
  std::cout << "Reading files..." << std::flush ;
  auto pointCloudInfo = read_pcfs_to_matrix(&m, &pcfs, (long)params.m_varSize, (long)params.m_offset);
  auto pointSize(pointCloudInfo.rows) ;
  double end(omp_get_wtime()) ;
  std::cout << "\t\t\t\t Done in " << end - start << "s \n" << std::endl ;

  std::cout << "File contains " << pointCloudInfo.rows << " rows and " << pointCloudInfo.columns << " columns. "
  << "Read data from columns " << (params.m_offset + 1) << " to " << (params.m_offset + params.m_varSize) << ".\n"
  << std::endl;

  // COMPUTING NORMALISED PROJECTION MATRIX
  start = omp_get_wtime();
  std::cout << "Computing projection matrix..." << std::flush ;
  MatrixXd pm = MatrixXd::Zero(timesSize, timesSize) ;
  pm = (1.0 / timesSize) * m.transpose() * m ;
  end = omp_get_wtime() ;
  std::cout << "\t\t\t Done in " << end - start << "s \n"
  << std::endl;

  // APPLY SPECTRAL POD FILTER IF DESIRED

  if (params.m_spodType > 0)
  {
    start = omp_get_wtime();
    std::cout << "Filtering projection matrix for SPOD..." << std::flush;

    int nfSize = 2 * params.m_spodWidth + 1;

    VectorXd g = VectorXd::Ones(nfSize);

    if (params.m_spodType == 2)
    {
      VectorXd gauss = VectorXd::LinSpaced(nfSize, -2.285, 2.285);
      g = exp(-square(gauss.array()));
    }

    g = g / g.sum();

    size_t idx = 0;

    MatrixXd spm = MatrixXd::Zero(timesSize, timesSize);
    MatrixXd pmExt = MatrixXd::Zero(timesSize * 3, timesSize * 3);

    pmExt = pm.replicate(3, 3).block(timesSize - params.m_spodWidth,
                                     timesSize - params.m_spodWidth,
                                     timesSize + 2 * params.m_spodWidth,
                                     timesSize + 2 * params.m_spodWidth);

    for (int i = 0; i < timesSize; i++)
    {
      for (int j = 0; j < timesSize; j++)
      {
        for (int k = -params.m_spodWidth; k < params.m_spodWidth + 1; k++)
        {
          spm(i, j) += g(idx) * pmExt(i + k + params.m_spodWidth, j + k + params.m_spodWidth);
          idx++;
        }
        idx = 0;
      }
    }

    pm = spm;

    end = omp_get_wtime();
    std::cout << "\t\t Done in " << end - start << "s \n"
    << std::endl;
  }

  // COMPUTING SORTED EIGENVALUES AND EIGENVECTORS

  start = omp_get_wtime();
  std::cout << "Computing eigenvalues and eigenvectors..." << std::flush;
  SelfAdjointEigenSolver<MatrixXd> eigensolver(pm);
  if (eigensolver.info() != Success)
    abort();

  VectorXd eigval = eigensolver.eigenvalues().reverse();
  MatrixXd eigvec = eigensolver.eigenvectors().rowwise().reverse();
  end = omp_get_wtime();
  std::cout << "\t Done in " << end - start << "s \n"
  << std::endl;

  // Adjust params.m_podSize according to the ric
  auto eigValSum(0.) ;
  for(auto i(0) ; i < eigval.size() ; ++i) {
    eigValSum += std::fabs(eigval(i)) ;
  }

  auto ric(0) ;
  auto podSizeWithRIC(0) ;
  while(ric < params.m_targetRic * eigValSum && podSizeWithRIC < eigval.size()) {
    ric += std::fabs(eigval(podSizeWithRIC)) ;
    ++podSizeWithRIC ;
  }
  params.m_podSize=std::min(podSizeWithRIC, params.m_podSize) ;
  std::cout << "With given RIC, pod size = " << params.m_podSize << std::endl;

  // COMPUTING POD MODES

  start = omp_get_wtime();
  std::cout << "Computing POD modes..." << std::flush;

  /* Define matrices to store scalar or vector values*/
  // For scalar or vector define one matrix
  const auto MVSIZE(pointSize*params.m_varSize) ;
  MatrixXd pod(MatrixXd::Zero(MVSIZE, params.m_podSize)) ;
  MatrixXd chronos(MatrixXd::Zero(params.m_podSize, timesSize)) ;

#pragma omp parallel
#pragma omp for
  for (size_t i = 0; i < params.m_podSize; i++)
  {
    for (size_t j = 0; j < timesSize; j++)
    {
      const auto factor(eigval(i) * timesSize) ;
      chronos(i, j) =  sqrt(factor) * eigvec(j, i) ;
      pod.col(i) += chronos(i, j) / factor * m.block(0, j, MVSIZE, 1);
    }
  }
  end = omp_get_wtime();
  std::cout << "\t\t\t\t Done in " << end - start << "s \n"
  << std::endl;

  // WRITING SORTED EIGENVALUES

  start = omp_get_wtime();
  std::cout << "Writing eigenvalues..." << std::flush;
  std::ofstream writeEigval(params.m_chronosDirName + "/eigenValues.bin", std::ios::binary);
  if (writeEigval.is_open()) {
    const auto size(eigval.size()) ;
    //writeEigval.write(reinterpret_cast<const char*>(&size), sizeof(size)) ;
    writeEigval.write(reinterpret_cast<const char*>(eigval.data()), size*sizeof(double)) ;
    writeEigval.close();
  }
  end = omp_get_wtime();
  std::cout << "\t\t\t\t Done in " << end - start << "s \n"
  << std::endl;

  // WRITING CHRONOS
  start = omp_get_wtime();
  std::cout << "Writing chronos..." << std::flush;
  std::ofstream writeChronos(params.m_chronosDirName + "/chronos.bin", std::ios::binary) ;
  if(writeChronos.is_open()) {
    //writeChronos.write(reinterpret_cast<const char*>(&params.m_podSize), sizeof(params.m_podSize)) ;
    //writeChronos.write(reinterpret_cast<const char*>(&timesSize), sizeof(timesSize)) ;
    //const auto chronosExtract(chronos.block(0, 0, params.m_podSize, timesSize)) ;
    writeChronos.write(reinterpret_cast<const char*>(chronos.data()), chronos.size() * sizeof(double)) ;
    writeChronos.close() ;
  }
  end = omp_get_wtime();
  std::cout << "\t\t\t\t Done in " << end - start << "s \n"
  << std::endl;

  // WRITING POD MODES
  start = omp_get_wtime();
  std::cout << "Writing POD modes..." << std::flush;
  std::ofstream writeMode(params.m_modeDirName + "/mode.bin", std::ios::binary);
  if(writeMode.is_open()) {
    const auto nbRows(pod.rows()) ;
    const auto nbCols(pod.cols()) ;
    //writeMode.write(reinterpret_cast<const char*>(&nbRows), sizeof(nbRows)) ;
    //writeMode.write(reinterpret_cast<const char*>(&nbCols), sizeof(nbCols)) ;
    writeMode.write(reinterpret_cast<const char*>(pod.data()), pod.size() * sizeof(double)) ;
    writeMode.close() ;
  }
  end = omp_get_wtime();
  std::cout << "\t\t\t\t Done in " << end - start << "s \n"
  << std::endl;
}

int main(int argc, const char *argv[])
{
  ez::ezOptionParser opt;

  opt.overview = "POD routine";
  opt.syntax = "Perform POD (Proper Orthogonal Decomposiztion) using [INPUTS] ...";
  opt.example = "Add example \n\n";
  opt.footer = "POD routine. \nThis program is free and without warranty.\n";

  opt.add(
      "",                            // Default.
      0,                             // Required?
      0,                             // Number of args expected.
      0,                             // Delimiter if expecting multiple args.
      "Display usage instructions.", // Help description.
      "-h"                          // Flag token.
      );

  opt.add(
      "",                             // Default.
      1,                              // Required?
      1,                              // Number of args expected.
      0,                              // Delimiter if expecting multiple args.
      "File with time snapshot list", // Help description.
      Parameters::m_timesFileNameOpt                          // Flag token.
      );

  opt.add(
      "",                                             // Default.
      1,                                              // Required?
      1,                                              // Number of args expected.
      0,                                              // Delimiter if expecting multiple args.
      "Point cloud file name (in time directories).", // Help description.
      Parameters::m_dataFileNameOpt                                         // Flag token.
      );

  opt.add(
      "",                                                                // Default.
      1,                                                                 // Required?
      1,                                                                 // Number of args expected.
      0,                                                                 // Delimiter if expecting multiple args.
      "Directory where the time directories reside as sub-directories.", // Help description.
      Parameters::m_inputDirNameOpt                                                              // Flag token.
      );

  opt.add(
      "",                   // Default.
      1,                    // Required?
      1,                    // Number of args expected.
      0,                    // Delimiter if expecting multiple args.
      "Chronos directory.", // Help description.
      Parameters::m_chronosDirNameOpt                 // Flag token.
      );

  opt.add(
      "",                 // Default.
      1,                  // Required?
      1,                  // Number of args expected.
      0,                  // Delimiter if expecting multiple args.
      "Modes directory.", // Help description.
      Parameters::m_modeDirNameOpt               // Flag token.
      );

  ez::ezOptionValidator *vS4 = new ez::ezOptionValidator("s4", "ge", "0");
  opt.add(
      "",                            // Default.
      1,                             // Required?
      1,                             // Number of args expected.
      0,                             // Delimiter if expecting multiple args.
      "Number of values per point.", // Help description.
      Parameters::m_varSizeOpt,                          // Flag token.
      vS4                            // Validate input
      );

  opt.add(
      "0",                                                 // Default.
      0,                                                   // Required?
      1,                                                   // Number of args expected.
      0,                                                   // Delimiter if expecting multiple args.
      "Point cloud file column offset (for reading data)", // Help description.
      Parameters::m_offsetOpt,                                               // Flag token.
      vS4                                                  // Validate input
      );

  opt.add(
      "",                          // Default.
      1,                           // Required?
      1,                           // Number of args expected.
      0,                           // Delimiter if expecting multiple args.
      "Number of modes to write.", // Help description.
      Parameters::m_podSizeOpt,                       // Flag token.
      vS4                          // Validate input
      );

  opt.add(
      "",                            // Default.
      1,                             // Required?
      1,                             // Number of args expected.
      0,                             // Delimiter if expecting multiple args.
      "Number of parallel threads.", // Help description.
      Parameters::m_threadsSizeOpt,                         // Flag token.
      vS4                            // Validate input
      );

  ez::ezOptionValidator* vD = new ez::ezOptionValidator("d", "gele", "0,1");

  opt.add(
      "0.9",                                             // Default.
      1,                                              // Required?
      1,                                              // Number of args expected.
      0,                                              // Delimiter if expecting multiple args.
      "RIC", // Help description.
      Parameters::m_targetRicOpt,                                         // Flag token.
      vD                            // Validate input
      );

  ez::ezOptionValidator *vS1 = new ez::ezOptionValidator("s1", "ge", "0");

  opt.add(
      "0", // 0 : No SPOD, 1 : SPOD with box filter, 2 : SPOD with Gauss filter
      0,
      1,
      0,
      "SPOD",
      Parameters::m_spodTypeOpt,
      vS1
      );

  opt.add(
      "5",
      0,
      1,
      0,
      "SPOD filter width",
      Parameters::m_spodWidthOpt,
      vS1
      );

  // Perform the actual parsing of the command line.
  opt.parse(argc, argv);

  if (opt.isSet("-h"))
  {
    Usage(opt);
    return 1;
  }

  // Perform validations of input parameters.
  //
  // Check if directories exist.
  std::array<std::string, 3> dirflags = {Parameters::m_inputDirNameOpt, Parameters::m_chronosDirNameOpt, Parameters::m_modeDirNameOpt};
  for (auto &dirflag : dirflags)
  {
    if (opt.isSet(dirflag))
    {
      std::string inputdir;
      struct stat info;
      opt.get(dirflag.c_str())->getString(inputdir);

      if (stat(inputdir.c_str(), &info) != 0)
      {
        std::cerr << "ERROR: " << inputdir << " does not exist.\n\n";
        return 1;
      }

      if (!(info.st_mode & S_IFDIR))
      { // S_ISDIR() doesn't exist on my windows
        std::cerr << "ERROR: " << inputdir << " is not a directory.\n\n";
        return 1;
      }
    }
  }

  std::vector<std::string> badOptions;
  int i;
  if (!opt.gotRequired(badOptions))
  {
    for (i = 0; i < badOptions.size(); ++i)
      std::cerr << "ERROR: Missing required option " << badOptions[i] << ".\n\n";

    Usage(opt);
    return 1;
  }

  if (!opt.gotExpected(badOptions))
  {
    for (i = 0; i < badOptions.size(); ++i)
      std::cerr << "ERROR: Got unexpected number of arguments for option " << badOptions[i] << ".\n\n";

    Usage(opt);
    return 1;
  }

  std::string firstArg;
  if (opt.firstArgs.size() > 0)
    firstArg = *opt.firstArgs[0];

  pod(opt);
  return 0;
}

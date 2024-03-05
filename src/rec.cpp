//
// Created by encheryg on 22/04/2022.
//

#include <iostream>
#include <Eigen/Dense>
#include "ezOptionParser.hpp"
#include <omp.h>
#include <sys/stat.h>

#include "utils.h"

void reconstruct(ez::ezOptionParser &opt) {
  std::cout << "Starting reconstruction routine " << std::endl ;

  Parameters params(opt) ;

  std::vector<std::string> t;

  t = read_timefile(params.m_timesFileName);
  long TSIZE = t.size();

  // READING INPUT FILES

  /* Build a string array with all file names to be read into the matrix */
  std::vector<std::string> pcfs;
  for (std::vector<std::string>::iterator it = t.begin(); it != t.end(); ++it)
  {
    std::string name_temp = params.m_inputDirName + "/" + *it + "/" + params.m_dataFileName ;
    pcfs.push_back(name_temp);
  }

  MatrixXd snapshots;
  double start(omp_get_wtime()) ;
  std::cout << "Reading snapshots files..." << std::flush;
  auto pointCloudInfo = read_pcfs_to_matrix(&snapshots, &pcfs, (long)params.m_varSize, (long)params.m_offset);
  auto REF_MSIZE = pointCloudInfo.rows;
  double end(omp_get_wtime());
  const auto snapsReadingTime(end - start) ;
  std::cout << "\t\t\t\t Done in " << snapsReadingTime << "s \n"
  << std::endl;

  std::cout << "File contains " << pointCloudInfo.rows << " rows and " << pointCloudInfo.columns << " columns. "
  << "Read data from columns " << (params.m_offset + 1) << " to " << (params.m_offset + params.m_varSize) << ".\n"
  << std::endl;

  // READING MODE FILES
  omp_set_num_threads(params.m_threadsSize);
  start = omp_get_wtime();
  const auto MVSIZE(pointCloudInfo.rows * params.m_varSize) ;
  std::cout << "Reading modes..." << std::flush;
  std::ifstream readMode(params.m_modeDirName + "/mode.bin", std::ios::binary) ;
  if (!readMode.is_open())
    throw "Could not open mode file" ;

  readMode.seekg(0, std::ios::end);
  const auto NSIZE(readMode.tellg() / MVSIZE) ;
  MatrixXd m(MVSIZE, NSIZE);
  readMode.seekg(0, std::ios::beg);
  readMode.read(reinterpret_cast<char*>(m.data()), m.size() * sizeof(double)) ;

  end = omp_get_wtime();
  const auto modesReadingTime(end - start) ;
  std::cout << "\t\t\t\t Done in " << modesReadingTime << "s \n"
  << std::endl;

  // COMPUTING BASES COEFFICIENTS
  start = omp_get_wtime();
  std::cout << "Computing coefficients..." << std::flush;
  MatrixXd c = MatrixXd::Zero(TSIZE, NSIZE);
#pragma omp parallel
#pragma omp for
  for (size_t l = 0; l < NSIZE; l++)
  {
    for (size_t k = 0; k < TSIZE; k++)
    {
      for (size_t j = 0; j < params.m_varSize ; j++){
        for(auto i = 0; i < REF_MSIZE; i++) {
          c(k, l) += m(i + REF_MSIZE * j, l) * snapshots(i + REF_MSIZE * j, k) ;
        }
      }
    }
  }
  end = omp_get_wtime();
  const auto coeffComputingTime(end - start) ;
  std::cout << "\t\t\t\t Done in " << coeffComputingTime << "s \n"
  << std::endl;

  // COMPUTE RECONSTRUCTED FIELDS
  start = omp_get_wtime();
  std::cout << "Computing reconstructed fields..." << std::flush;
  MatrixXd rec = MatrixXd::Zero(MVSIZE, TSIZE);
#pragma omp parallel
#pragma omp for
  for (size_t i = 0; i < TSIZE; i++)
  {
    for (size_t j = 0; j < NSIZE; j++)
    {
      rec.col(i) = rec.col(i) + c(i, j) * m.block(0, j, MVSIZE, 1);
    }
  }
  end = omp_get_wtime();
  const auto recComputingTime(end - start) ;
  std::cout << "\t\t Done in " << recComputingTime << "s \n" << std::endl;

  // WRITE RECONSTRUCTED FIELDS
  start = omp_get_wtime();
  std::cout << "Writing reconstructed fields..." << std::flush;
  std::ofstream writeField(params.m_recDirName + "/reconstruction.bin", std::ios::binary);
  if(writeField.is_open()) {
    writeField.write(reinterpret_cast<const char*>(rec.data()), rec.size() * sizeof(double)) ;
    writeField.close() ;
  }
  end = omp_get_wtime();
  const auto resWritingTime(end - start) ;
  std::cout << "\t\t\t Done in " << resWritingTime << "s \n"
  << std::endl;

  const auto globalTime(snapsReadingTime+modesReadingTime+coeffComputingTime+recComputingTime+resWritingTime) ;
  std::cout << "Everything done in " << globalTime << "s \n" << std::endl;
}

int main(int argc, const char *argv[])
{
  ez::ezOptionParser opt;

  opt.overview = "Reconstruction routine";
  opt.syntax = "Perform a reconstruction based on POD modes using [INPUTS] ...";
  opt.example = "Add example \n\n";
  opt.footer = "\nThis program is free and without warranty.\n";

  opt.add(
      "",                            // Default.
      0,                             // Required?
      0,                             // Number of args expected.
      0,                             // Delimiter if expecting multiple args.
      "Display usage instructions.", // Help description.
      "-h"                          // Flag token.
      );

  ez::ezOptionValidator *vS4 = new ez::ezOptionValidator("s4", "ge", "0");
  opt.add(
      "",                            // Default.
      1,                             // Required?
      1,                             // Number of args expected.
      0,                             // Delimiter if expecting multiple args.
      "Number of values per point.", // Help description.
      Parameters::m_varSizeOpt,                          // Flag token.
      vS4                            //Validate input
      );

  opt.add(
      "0",                                                 // Default.
      0,                                                   // Required?
      1,                                                   // Number of args expected.
      0,                                                   // Delimiter if expecting multiple args.
      "Point cloud file column offset (for reading data)", // Help description.
      "-co",                                               // Flag token.
      Parameters::m_offsetOpt,
      vS4                           // Validate input
      );

  opt.add(
      "",                            // Default.
      1,                             // Required?
      1,                             // Number of args expected.
      0,                             // Delimiter if expecting multiple args.
      "Number of parallel threads.", // Help description.
      Parameters::m_threadsSizeOpt,                         // Flag token.
      vS4                            //Validate input
      );

  opt.add(
      "",                                                                // Default.
      1,                                                                 // Required?
      1,                                                                 // Number of args expected.
      0,                                                                 // Delimiter if expecting multiple args.
      "Directory where the time directories reside as sub-directories.", // Help description.
      Parameters::m_inputDirNameOpt // Flag token.
      );

  opt.add(
      "",                 // Default.
      1,                  // Required?
      1,                  // Number of args expected.
      0,                  // Delimiter if expecting multiple args.
      "Reconstruction output directory.", // Help description.
      Parameters::m_recDirNameOpt               // Flag token.
      );

  opt.add(
      "",                 // Default.
      1,                  // Required?
      1,                  // Number of args expected.
      0,                  // Delimiter if expecting multiple args.
      "Modes directory.", // Help description.
      Parameters::m_modeDirNameOpt               // Flag token.
      );

  opt.add(
      "",                             // Default.
      1,                              // Required?
      1,                              // Number of args expected.
      0,                              // Delimiter if expecting multiple args.
      "File with time snapshot list", // Help description.
      Parameters::m_timesFileNameOpt                           // Flag token.
      );

  opt.add(
      "",                                             // Default.
      1,                                              // Required?
      1,                                              // Number of args expected.
      0,                                              // Delimiter if expecting multiple args.
      "Point cloud file name (in time directories).", // Help description.
      Parameters::m_dataFileNameOpt// Flag token.
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
  std::array<std::string, 3> dirflags = {Parameters::m_inputDirNameOpt, Parameters::m_modeDirNameOpt,
                                         Parameters::m_recDirNameOpt};
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

  reconstruct(opt) ;

  return 0;
}
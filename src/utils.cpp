//
// Created by encheryg on 02/05/2022.
//

#include "utils.h"

std::vector<std::string> read_timefile(const std::string tfile)
{
  /* Open and read file with list of times */
  std::ifstream timefile(tfile);
  std::vector<std::string> t_entr;

  if (timefile.is_open())
  {
    /* Count number of words in file. This should correspond to the
    number of snapshots to use. */
    std::string line;

    while (getline(timefile, line))
    {
      t_entr.push_back(line);
    }
  }
  timefile.close();

  return t_entr;
}

/*
Function for splitting strings into a vector datatype.
*/
std::vector<std::string> split_string(const std::string &s, char delimiter)
{
  std::vector<std::string> tokens;
  std::string token;
  std::istringstream tokenStream(s);
  while (std::getline(tokenStream, token, delimiter))
  {
    tokens.push_back(token);
  }
  return tokens;
}

/*
Parse the point cloud files and populate matrix with data.
*/
pointCloudFileInfo read_pcfs_to_matrix(MatrixXd *m,
                                       const std::vector<std::string> *fvec,
                                       const long no_cols,
                                       const long offset)
{
  pointCloudFileInfo pointCloudRefFileInfo;
  bool verbose = false;

  if (verbose)
  {
    std::ostream_iterator<std::string> out_it(std::cout, "\n");
    std::copy(fvec->begin(), fvec->end(), out_it);
  }

  /* Establish a reference number of points for checking the problem size.
  The size is determined from the point cloud file in the first time directory. */
  pointCloudRefFileInfo.rows = 0;
  std::string ref_fname = *(fvec->begin());

  std::ifstream ref_file(ref_fname);
  if (ref_file.is_open())
  {
    std::string unused_line;

    /* Get a line and count the number of columns by splitting at spaces ' '.
    Increment the row counter. */
    getline(ref_file, unused_line);
    pointCloudRefFileInfo.rows++;
    auto tokens = split_string(unused_line, ' ');
    pointCloudRefFileInfo.columns = tokens.size();

    /* Count the rest of the rows. */
    while (getline(ref_file, unused_line))
      pointCloudRefFileInfo.rows++;
  }
  ref_file.close();

  if (verbose)
  {
    std::cout << "From file " << ref_fname << std::endl;
    std::cout << "found " << pointCloudRefFileInfo.rows << " points" << std::endl;
  }

  /* Number of time samples to consider is determined from the length of the
  vector of files */
  long TSIZE = fvec->size();

  /* Define matrix to store the file content */
  *m = MatrixXd::Zero(pointCloudRefFileInfo.rows * no_cols, TSIZE);
  double dummy = 0;
#pragma omp parallel
#pragma omp for
  for (size_t snapshot = 0; snapshot < TSIZE; snapshot++)
  {
    std::ifstream file((*fvec)[snapshot]);

    if (file.is_open())
    {
      for (size_t row_idx = 0; row_idx < pointCloudRefFileInfo.rows; row_idx++)
      {
        for (size_t col_no = 0, j = 0; col_no < (no_cols + offset); col_no++)
        {
          if (col_no < offset)
            file >> dummy;
          else
          {
            file >> (*m)(row_idx + pointCloudRefFileInfo.rows * j, snapshot);
            j++;
          }
        }
      }
      file.close();
    }
    else
    {
      std::cerr << "Unable to open file " << (*fvec)[snapshot] << std::endl;
    }
  }

  return pointCloudRefFileInfo;
}

void Usage(ez::ezOptionParser &opt)
{
  std::string usage;
  opt.getUsage(usage);
  std::cout << usage;
};

const char* Parameters::m_varSizeOpt = "-v";
const char* Parameters::m_offsetOpt = "-co" ;
const char* Parameters::m_podSizeOpt = "-nm" ;
const char* Parameters::m_threadsSizeOpt = "-np" ;
const char* Parameters::m_inputDirNameOpt = "-i" ;
const char* Parameters::m_chronosDirNameOpt = "-c" ;
const char* Parameters::m_modeDirNameOpt = "-m";
const char* Parameters::m_timesFileNameOpt = "-tf" ;
const char* Parameters::m_dataFileNameOpt = "-pcfn" ;
const char* Parameters::m_targetRicOpt = "-ric";
const char* Parameters::m_spodTypeOpt = "-spod-type" ;
const char* Parameters::m_spodWidthOpt = "-spod-width" ;
const char* Parameters::m_recDirNameOpt = "-r" ;



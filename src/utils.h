//
// Created by encheryg on 08/04/2022.
//

#ifndef POD_UTILS_H
#define POD_UTILS_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <Eigen/Dense>
#include <iterator>

#include "ezOptionParser.hpp"
#include <Eigen/Dense>

/*
Function for reading the input file containing all the itme entries to use.
*/
std::vector<std::string> read_timefile(const std::string tfile) ;

/*
Function for splitting strings into a vector datatype.
*/
std::vector<std::string> split_string(const std::string &s, char delimiter) ;

/*
Struct to hold information about the point cloud files.
*/
struct pointCloudFileInfo {
  long rows;
  long columns;
};

using namespace Eigen;

/*
Parse the point cloud files and populate matrix with data.
*/
pointCloudFileInfo read_pcfs_to_matrix(MatrixXd *m,
                                       const std::vector<std::string> *fvec,
                                       const long no_cols,
                                       const long offset) ;

void Usage(ez::ezOptionParser &opt) ;

struct Parameters {
  Parameters(ez::ezOptionParser &opt) :
  m_varSize(0),
  m_offset(0),
  m_podSize(0),
  m_threadsSize(0),
  m_inputDirName(""),
  m_chronosDirName(""),
  m_modeDirName(""),
  m_timesFileName(""),
  m_dataFileName(""),
  m_recDirName(""),
  m_targetRic(0.),
  m_spodType(0),
  m_spodWidth(0) {
    if(opt.isSet(m_varSizeOpt))
      opt.get(m_varSizeOpt) -> getInt(m_varSize) ;

    if(opt.isSet(m_offsetOpt))
      opt.get(m_offsetOpt) -> getInt(m_offset) ;

    if(opt.isSet(m_podSizeOpt))
      opt.get(m_podSizeOpt) -> getInt(m_podSize) ;

    if(opt.isSet(m_threadsSizeOpt))
      opt.get(m_threadsSizeOpt) -> getInt(m_threadsSize) ;

    if(opt.isSet(m_inputDirNameOpt))
      opt.get(m_inputDirNameOpt) -> getString(m_inputDirName) ;

    if(opt.isSet(m_chronosDirNameOpt))
      opt.get(m_chronosDirNameOpt) -> getString(m_chronosDirName) ;

    if(opt.isSet(m_modeDirNameOpt))
      opt.get(m_modeDirNameOpt) -> getString(m_modeDirName) ;

    if(opt.isSet(m_timesFileNameOpt))
      opt.get(m_timesFileNameOpt) -> getString(m_timesFileName) ;

    if(opt.isSet(m_dataFileNameOpt))
      opt.get(m_dataFileNameOpt) -> getString(m_dataFileName) ;

    if(opt.isSet(m_recDirNameOpt))
      opt.get(m_recDirNameOpt) -> getString(m_recDirName) ;

    if(opt.isSet(m_targetRicOpt))
      opt.get(m_targetRicOpt) -> getDouble(m_targetRic) ;

    if(opt.isSet(m_spodTypeOpt))
      opt.get(m_spodTypeOpt) -> getInt(m_spodType) ;

    if(opt.isSet(m_spodWidthOpt))
      opt.get(m_spodWidthOpt) -> getInt(m_spodWidth) ;
  }

  int m_varSize ;
  int m_offset ;
  int m_podSize ;
  int m_threadsSize ;
  std::string m_inputDirName ;
  std::string m_chronosDirName ;
  std::string m_modeDirName ;
  std::string m_timesFileName ;
  std::string m_dataFileName ;
  std::string m_recDirName ;
  double m_targetRic ;
  int m_spodType ;
  int m_spodWidth ;

  static const char* m_varSizeOpt ;
  static const char* m_offsetOpt ;
  static const char* m_podSizeOpt ;
  static const char* m_threadsSizeOpt ;
  static const char* m_inputDirNameOpt ;
  static const char* m_chronosDirNameOpt ;
  static const char* m_modeDirNameOpt ;
  static const char* m_timesFileNameOpt ;
  static const char* m_dataFileNameOpt ;
  static const char* m_targetRicOpt ;
  static const char* m_spodTypeOpt ;
  static const char* m_spodWidthOpt ;
  static const char* m_recDirNameOpt ;
} ;

#endif //POD_UTILS_H

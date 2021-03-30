/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2020 Triad National Security, LLC
 * All rights reserved.
 * --------------------------------------------------------------------------~*/

/*~--------------------------------------------------------------------------~*
 *
 * /@@@@@@@@  @@           @@@@@@   @@@@@@@@ @@@@@@@  @@      @@
 * /@@/////  /@@          @@////@@ @@////// /@@////@@/@@     /@@
 * /@@       /@@  @@@@@  @@    // /@@       /@@   /@@/@@     /@@
 * /@@@@@@@  /@@ @@///@@/@@       /@@@@@@@@@/@@@@@@@ /@@@@@@@@@@
 * /@@////   /@@/@@@@@@@/@@       ////////@@/@@////  /@@//////@@
 * /@@       /@@/@@//// //@@    @@       /@@/@@      /@@     /@@
 * /@@       @@@//@@@@@@ //@@@@@@  @@@@@@@@ /@@      /@@     /@@
 * //       ///  //////   //////  ////////  //       //      //
 *
 *~--------------------------------------------------------------------------~*/

/**
 * @file influx.h
 * @brief Routines for flux injection, specifically for kilonova ejecta
 *
 * This file contains a set of tools to read and interpolate input flux from
 * other codes, particularly from nubhlight.
 *
 * @author Oleg Korobkin, Hyun Lim
 * @date March 2021
 */
#pragma once

#include <iostream>
#include <cmath>
#include <vector>
#include <boost/algorithm/string.hpp>
#include <fstream>
#include <cstdio>
#include "eos.h"
#include "body.h"
#include "params.h"

#include <glob.h>

namespace influx {

static std::vector<std::string> input_filenames;
static std::vector<double> grid_times;
static std::vector<double> grid_theta;
static std::vector<double> grid_phi;

/**
* @brief   Returns a list of file names matching the pattern
*
* From stackoverflow:
* https://stackoverflow.com/questions/8401777/simple-glob-in-c-on-unix-system
*/
void
glob_input_filenames(const std::string& pattern) {
  using namespace std;

  // glob result
  glob_t glob_result;

  // wipe out before use
  memset(&glob_result, 0, sizeof(glob_t));

  // do the glob
  int retval = glob(pattern.c_str(), GLOB_TILDE, NULL, &glob_result);
  if (retval != 0) {
    globfree(&glob_result);
    cerr << "ERROR: unable to expand pattern " << pattern << "!" << endl;
    exit(1);
  }

  // collect all the filenames
  input_filenames.empty();
  for (size_t i = 0; i < glob_result.gl_pathc; ++i) {
    input_filenames.push_back(string(glob_result.gl_pathv[i]));
  }

}

/**
* @brief   Reads a time stamp from a single flux file
*/
double
read_time_stamp(const std::string & filename) {
  using namespace std;
  ifstream infile;
  string line;

  // attempt to open the file
  infile.open(filename.c_str());
  if(!infile) {
    cerr << "ERROR: Unable to open file '" << filename << "'" << endl;
    exit(1);
  }

  // read 4th line
  for(int ln = 1; ln<=4; ++ln) {
    std::getline(infile, line);
  }
  infile.close();

  size_t eqpos = line.find("=");
  if (eqpos == string::npos) {
    cerr << "ERROR: cannot find time stamp in 4th line of '"
         << filename << "'" << endl;
    exit(1);
  }
  istringstream iss(line.substr(eqpos + 1));
  double t = 0.0;
  iss >> t;
  return t;

}

/**
* @brief   Reads all time stamps
*
* Uses: grid_times
*/
void
read_time_stamps() {
  using namespace std;
  grid_times.empty();

  double t = 0.0;
  for (size_t i = 0; i < input_filenames.size(); ++i) {
    string fname = input_filenames[i];
    double tm = read_time_stamp(fname);
    grid_times.push_back(tm);
  }
}

/**
* @brief   Reads spherical grid from a file
*
* Expects a file in ASCII format (gnuplot 2D data format),
* with several lines of comments/ metadata in the header
* (starting with '#' symbol). The data is organized in several
* columns (at least three), separated by spaces, and in blocks,
* separateb by an empty line:
*
* -- >>> flux input file example >>> ---------------------------
* # commment comment
* # metadata
* # 1:theta 2:phi      3:rho[g/cm3]   4:v_r[c]
* 0.000     0.000      0.000e+00      0.000e+00
* 0.000     0.100      0.000e+00      0.000e+00
* 0.000     0.200      0.000e+00      0.000e+00
* 0.000     0.300      0.000e+00      0.000e+00
*
* 1.000     0.000      0.000e+00      0.000e+00
* 1.000     0.100      0.000e+00      0.000e+00
* 1.000     0.200      0.000e+00      0.000e+00
* 1.000     0.300      0.000e+00      0.000e+00
*
* 2.000     0.000      0.000e+00      0.000e+00
* 2.000     0.100      0.000e+00      0.000e+00
* 2.000     0.200      0.000e+00      0.000e+00
* 2.000     0.300      0.000e+00      0.000e+00
* -- <<< -------------------------------------------------------
*
* Modifies namespace variables:
*  - grid_theta : number of grid points in theta direction
*  - grid_phi   : number of grid points in the phi-direction
*
*/
void
read_spherical_grid(const std::string & filename) {
  using namespace std;
  ifstream infile;
  string line;

  // attempt to open the file
  infile.open(filename.c_str());
  if(!infile) {
    cerr << "ERROR: Unable to open file '" << filename << "'" << endl;
    exit(1);
  }

  // skip the header or blank lines
  double theta = -1.0;
  double phi = -1.0;
  for(int ln = 1; std::getline(infile, line); ++ln) {

    // skip comments (lines starting with '#' at any position)
    // and blank lines
    bool is_blank = true, is_comment = false;
    for(size_t i = 0; i < line.length(); i++) {
      char c = line[i];
      is_comment = (c == '#');
      is_blank = (c == ' ' || c == '\t');
      if(is_comment or not is_blank)
        break;
    }
    if(is_comment or is_blank)
      continue;

    // read two values from line, add to phi/theta grids if new
    istringstream iss(line);
    double theta_next, phi_next;
    iss >> theta_next;
    iss >> phi_next;
    if (theta_next > theta) {
      theta = theta_next;
      grid_theta.push_back(theta);
    }
    if (phi_next > phi) {
      phi = phi_next;
      grid_phi.push_back(phi);
    }

  }
  infile.close();

}

/**
* @brief   Initial scan of the flux data and parameter setup
*
*/
void
init() {
  using namespace std;
  using namespace param;
  string line;

  glob_input_filenames(input_flux_files);
  sort(input_filenames.begin(), input_filenames.end());
  read_time_stamps();
  read_spherical_grid(input_filenames[0]);

  for (auto t : grid_times) { cout << t << " "; }
  cout << endl;
  for (auto t : grid_theta) { cout << t << " "; }
  cout << endl;
  for (auto t : grid_phi) { cout << t << " "; }
  cout << endl;
  exit(0);
}



} // namespace influx


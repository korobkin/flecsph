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

static std::vector<double> grid_times;  // 1D grid of all timesteps
static std::vector<double> grid_theta;  // 1D grid of theta: varies with t!
static std::vector<double> grid_phi;    // 1D grid of phi
static size_t INFLX_NT = 0;             // total number of timesteps
static size_t INFLX_NT_WINDOW = 0;      // timestep window
static size_t INFLX_NTHETA = 0;         // number of grid points in theta-direction
static size_t INFLX_NPHI = 0;           // number of grid points in phi-direction
static double extraction_radius = 0.0;  // extraction radius (read from the files)

#define IND3(IT,JTH,KPHI) ((KPHI)+INFLX_NPHI*((JTH)+INFLX_NTHETA*(IT)))
static std::vector<double> grid3d_theta;
static std::vector<double> grid3d_rho;
static std::vector<double> grid3d_vr;
static std::vector<double> grid3d_vth;
static std::vector<double> grid3d_vphi;
static std::vector<double> grid3d_uint;
static std::vector<double> grid3d_pres;
static std::vector<double> grid3d_temp;
static std::vector<double> grid3d_ye;

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
* 
* The time stamp is in the 4th line, looking e.g. like this:
* # Time [s]: t= 0.15889199090321768
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
    if (tm < t) {
      cerr << "ERROR: flux files not in order for files "
           << "'" << input_filenames[i-1] << "', "
           << "'" << fname << endl;
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
    grid_times.push_back(tm);
    t= tm;
  }
}

/**
* @brief   Reads an extraction radius from a single flux file
* 
* The extraction radius is in the 5th line:
* # Extraction radius [cm]: R= 349828765.9436025 
*/
double
read_extraction_radius(const std::string & filename) {
  using namespace std;
  ifstream infile;
  string line;

  // attempt to open the file
  infile.open(filename.c_str());
  if(!infile) {
    cerr << "ERROR: Unable to open file '" << filename << "'" << endl;
    exit(1);
  }

  // read 5th line
  for(int ln = 1; ln<=5; ++ln) {
    std::getline(infile, line);
  }
  infile.close();

  size_t eqpos = line.find("R=");
  if (eqpos == string::npos) {
    cerr << "ERROR: cannot find time stamp in 5th line of '"
         << filename << "'" << endl;
    exit(1);
  }
  istringstream iss(line.substr(eqpos + 2));
  double R_ext = 0.0;
  iss >> R_ext;
  return R_ext;

} // read_extraction_radius

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

} // read_spherical_grid()

/**
* @brief   Reads the density and velocity from a single flux file
*
* See flux input file format in 'read_spherical_grid' function
* Parameters:
*  - filename:   file to read
*  - it:         time index in data arrays where to store the data
*  - all_vars:   read all variables (otherwise, only rho, v and theta)
*/
void
read_single_snap(const std::string & filename, 
    const int it, const bool all_vars) {
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
  int n = 0, ith, jphi;
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

    ith  = n / INFLX_NPHI;
    jphi = n % INFLX_NPHI;

    // read theta, phi, rho and vr
    istringstream iss(line);
    double theta, phi, rho, vr;
    iss >> theta >> phi >> rho >> vr;

    grid3d_theta[IND3(it,ith,jphi)] = theta;
    grid3d_rho[IND3(it,ith,jphi)] = rho;
    grid3d_vr[IND3(it,ith,jphi)] = vr*C_LIGHT_CGS;

    if (all_vars) {
      double vth, vphi, ye, temp, uint, pres;
      iss >> vth >> vphi >> uint >> pres >> temp >> ye;
      grid3d_vth[IND3(it,ith,jphi)] = vth*C_LIGHT_CGS;
      grid3d_vphi[IND3(it,ith,jphi)] = vphi*C_LIGHT_CGS;
      grid3d_uint[IND3(it,ith,jphi)] = uint;
      grid3d_pres[IND3(it,ith,jphi)] = pres;
      grid3d_temp[IND3(it,ith,jphi)] = temp;
      grid3d_ye[IND3(it,ith,jphi)] = ye;
    }
    ++n;

  }
  infile.close();

} // read_single_snap

/**
* @brief   Reads all time stamps
*
* Uses: grid_times
*/
double
compute_total_mass() {
  using namespace std;
  double mass = 0.0;
  const double dphi = 2.0*M_PI/(double)INFLX_NPHI;
  std::vector<double> dth;

  for (int it=1; it<INFLX_NT; ++it) {
    read_single_snap(input_filenames[it], it, true);
    double dt = grid_times[it] - grid_times[it-1];

    // compute angular differentials (different for every timestep)
    // TODO: remove the third index
    dth.empty();
    for (int ith = 0; ith < INFLX_NTHETA-1; ++ith) {
      dth.push_back(grid3d_theta[IND3(it,ith+1,0)]
                  - grid3d_theta[IND3(it,ith,0)]);
    }
    dth.push_back(grid3d_theta[IND3(it,0,0)] 
                - grid3d_theta[IND3(it,INFLX_NTHETA-1,0)]
                + M_PI);
    for (int ith = 0; ith < INFLX_NTHETA; ++ith) {
      double dm = 0.0;
      for (int jphi = 0; jphi < INFLX_NPHI; ++jphi) {
        double vr = grid3d_vr[IND3(it,ith,jphi)];
        if (vr > 0) 
          dm += grid3d_rho[IND3(it,ith,jphi)]*vr;
      }
      mass += dm*dt*dth[ith]*sin(grid3d_theta[IND3(it,ith,0)]);
    }
  }
  mass *= dphi*extraction_radius*extraction_radius;
  return mass;

} // compute_total_mass

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

  // consistency check: see that all extraction radii are identical
  extraction_radius = read_extraction_radius(input_filenames[0]);
  for(auto fname : input_filenames) {
    double R_ext = read_extraction_radius(fname);
    mpi_assert(R_ext == extraction_radius);
  }

  // set grid size
  read_spherical_grid(input_filenames[0]);
  INFLX_NT = grid_times.size();
  INFLX_NPHI = grid_phi.size();
  INFLX_NTHETA = grid_theta.size();
  INFLX_NT_WINDOW = 1;

  // resize data arrays
  const size_t N3D = INFLX_NT*INFLX_NPHI*INFLX_NTHETA;
  grid3d_theta.resize(N3D);
  grid3d_rho.resize(N3D);
  grid3d_vr.resize(N3D);
  grid3d_vth.resize(N3D);
  grid3d_vphi.resize(N3D);
  grid3d_uint.resize(N3D);
  grid3d_pres.resize(N3D);
  grid3d_temp.resize(N3D);
  grid3d_ye.resize(N3D);
  double mass = compute_total_mass();

  cout << "total mass: " << (mass/M_SUN_CGS) << endl;
  exit(0);
}

#undef IND3

} // namespace influx


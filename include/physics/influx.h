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
#include <hdf5.h>
#include "h5aux.h"

#include <glob.h>

namespace influx {

static std::vector<std::string> input_filenames;

// handy index macro (defined only inside this namespace)
#define IND2(IT,JTH) ((JTH)+INFLX_NTHETA*(IT))
#define IND3(IT,JTH,KPHI) ((KPHI)+INFLX_NPHI*((JTH)+INFLX_NTHETA*(IT)))

template <typename T>
struct grid_data_point_u {
  // fields with units as they appear in flux data files
  T   rho,
       vr,
      vth,
     vphi,
     uint,
     pres,
     temp,
       ye;


  //! Default constructor.
  grid_data_point_u() : rho(0.), vr(0.), vth(0.), vphi(0.),
      uint(0.), pres(0.), temp(0.), ye(0.) {}

  //! Default copy constructor.
  grid_data_point_u(grid_data_point_u const &) = default;

  //! Assignment operator.
  grid_data_point_u & operator=(grid_data_point_u const & rhs) {
    if(this != &rhs) {
      rho = rhs.rho;
      vr = rhs.vr;
      vth = rhs.vth;
      vphi = rhs.vphi;
      uint = rhs.uint;
      pres = rhs.pres;
      temp = rhs.temp;
      ye = rhs.ye;
    } // if

    return *this;
  } // operator =

  //--------------------------------------------------------------------------//
  // Macro to avoid code replication.
  //--------------------------------------------------------------------------//

#define def_operator(op)                                                       \
  grid_data_point_u & operator op(grid_data_point_u const & rhs) {             \
    rho op rhs.rho;                                                            \
    vr op rhs.vr;                                                              \
    vth op rhs.vth;                                                            \
    vphi op rhs.vphi;                                                          \
    uint op rhs.uint;                                                          \
    pres op rhs.pres;                                                          \
    temp op rhs.temp;                                                          \
    ye op rhs.ye;                                                              \
    return *this;                                                              \
  }

  def_operator(+=)
  def_operator(-=)
  def_operator(*=)
  def_operator(/=)

#define def_operator_type(op)                                                  \
  grid_data_point_u & operator op(T val) {                                     \
    rho op val;                                                                \
    vr op val;                                                                 \
    vth op val;                                                                \
    vphi op val;                                                               \
    uint op val;                                                               \
    pres op val;                                                               \
    temp op val;                                                               \
    ye op val;                                                                 \
    return *this;                                                              \
  }
  def_operator_type(+=)
  def_operator_type(-=)
  def_operator_type(*=)
  def_operator_type(/=)

};

/*!
  \function      operator+(grid_data_point_u, grid_data_point_u)
  \brief         Addition operator between two grid points

  \tparam T      Data type
 */
template<typename T>
grid_data_point_u<T>
operator+(const grid_data_point_u<T> & a, const grid_data_point_u<T> & b) {
  grid_data_point_u<T> tmp(a);
  tmp += b;
  return tmp;
} // operator +

/*!
  \function      operator-(grid_data_point_u, grid_data_point_u)
  \brief         Subtraction operator between two grid points

  \tparam T      Data type
 */
template<typename T>
grid_data_point_u<T>
operator-(const grid_data_point_u<T> & a, const grid_data_point_u<T> & b) {
  grid_data_point_u<T> tmp(a);
  tmp -= b;
  return tmp;
} // operator -

/*!
  \function      operator*(grid_data_point_u, T)
  \brief         Multiplication by a scalar

  \tparam T      Data type
 */
template<typename T>
grid_data_point_u<T>
operator*(const grid_data_point_u<T> & a, const T & b) {
  grid_data_point_u<T> tmp(a);
  tmp *= b;
  return tmp;
} // operator *

/*!
  \function      operator*(grid_data_point_u, T)
  \brief         Multiplication by a scalar

  \tparam T      Data type
 */
template<typename T>
grid_data_point_u<T>
operator*(const T & b, const grid_data_point_u<T> & a) {
  grid_data_point_u<T> tmp(a);
  tmp *= b;
  return tmp;
} // operator *

using grid_data_point_t = grid_data_point_u<double>;


static std::vector<double> grid_times;  // 1D grid of all timesteps
static std::vector<double> grid_theta;  // 1D grid of theta: varies with t!
static std::vector<double> grid_phi;    // 1D grid of phi
static size_t INFLX_NT = 0;             // total number of timesteps
static size_t INFLX_NTHETA = 0;         // number of grid points in theta-direction
static size_t INFLX_NPHI = 0;           // number of grid points in phi-direction
static double extraction_radius = 0.0;  // [cm] extraction radius (read from the files)
static double total_ejecta_mass = 0.0;  // [Msun] total mass

static std::vector<grid_data_point_t> grid3d_data;
static std::vector<double> grid2d_theta, grid2d_dth;
static hsize_t grid3d_dims[3];

// The first quantity records partially summed mass for each cell of the grid
// and the second one is a 1D array of partially summed mass up time step
static std::vector<double> grid3d_cumulative_mass, grid1d_cumulative_mass;

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
* Returns: a pair [Ntheta, Nphi] of dimensions
*
*/
std::pair<size_t, size_t>
read_spherical_grid_dimensions(const std::string & filename) {
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

  size_t N_theta, N_phi;
  N_theta = N_phi = 0;
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
      N_theta++;
    }
    if (phi_next > phi) {
      phi = phi_next;
      N_phi++;
    }

  }
  infile.close();

  return {N_theta, N_phi};

} // read_spherical_grid_dimensions()

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
    auto & gp = grid3d_data[IND3(it,ith,jphi)]; // gp - grid point
    istringstream iss(line);
    double phi; // dummy
    iss >> grid2d_theta[IND2(it,ith)] >> phi >> gp.rho >> gp.vr;
    gp.vr *= C_LIGHT_CGS;

    if (all_vars) {
      iss >> gp.vth >> gp.vphi >> gp.uint >> gp.pres >> gp.temp >> gp.ye;
      gp.vth  *= C_LIGHT_CGS;
      gp.vphi *= C_LIGHT_CGS;
    }
    ++n;

  }
  infile.close();

  // shift the theta array to span from 0 to PI-dth:
  size_t ij = IND2(it,INFLX_NTHETA-1);
  for (ith=INFLX_NTHETA-1; ith>0; --ith, --ij) {
    grid2d_theta[ij] *= 0.5;
    grid2d_theta[ij] += 0.5*grid2d_theta[ij-1];
  }
  grid2d_theta[IND2(it,0)] = 0.0;

  // fill array dtheta:
  ij = IND2(it,1);
  for (ith=0; ith<INFLX_NTHETA-1; ++ith,++ij) {
    grid2d_dth[ij] = grid2d_theta[ij+1] - grid2d_theta[ij];
  }
  ij = IND2(it,INFLX_NTHETA-1);
  grid2d_dth[ij] = M_PI - grid2d_theta[ij];

} // read_single_snap

/**
* @brief   Integrates the ejected mass over all timesteps
*/
double
compute_total_mass() {
  using namespace std;
  double mass = 0.0;
  const double dphi = 2.0*M_PI/(double)INFLX_NPHI;
  const double d3 = extraction_radius*dphi;

  for (int it=1; it<INFLX_NT; ++it) {
    read_single_snap(input_filenames[it], it, true);
    double dt = grid_times[it] - grid_times[it-1];
    for (int ith = 0; ith < INFLX_NTHETA; ++ith) {
      double d2 = grid2d_dth[IND2(it,ith)]*extraction_radius;
      double sin_th = sin(grid2d_theta[IND2(it,ith)]);
      for (int jphi = 0; jphi < INFLX_NPHI; ++jphi) {
        auto gp = grid3d_data[IND3(it,ith,jphi)];
        double d1 = (gp.vr > 0) ? gp.vr*dt : 0.0;
        mass += gp.rho*d1*d2*d3*sin_th/M_SUN_CGS;
        grid3d_cumulative_mass[IND3(it-1,ith,jphi)] = mass;
      }
    }
    grid1d_cumulative_mass[it-1] = mass;
  }
  grid1d_cumulative_mass[INFLX_NT - 1] = mass;
  return mass;

} // compute_total_mass

/**
* @brief   Integrates the ejected mass over all timesteps
*
* VS:  template parameter must be a 'vector space', with
*      the addition and multiplication by scalar operations
*/
grid_data_point_t
linear_interpolator(const double tm, const double theta,
    const double phi) {
  const size_t it  = interp::get_index(tm, grid_times);
  double * theta_it = grid2d_theta.data() + it*INFLX_NTHETA;
  double * dth_it = grid2d_dth.data() + it*INFLX_NTHETA;
  const size_t jth = interp::get_index(theta, theta_it, INFLX_NTHETA);
  const double dphi = 2.*M_PI/(double)INFLX_NPHI;
  const size_t kphi = int(phi/dphi);

  grid_data_point_t retval;
  if (0 <= it && it < INFLX_NT-1) { // otherwise, return 0

    const size_t
      it1 = it + 1,
      jth1 = jth + 1,
      kphi1 = kphi + 1;

    const double
      theta_i = theta_it[jth],
      phi_i = kphi*dphi,
      tm_i = grid_times[it];

    const double
      f1 = (tm - tm_i)/(grid_times[it1] - tm_i),
      f0 = 1. - f1;

    const double
      g1 = (theta - theta_i)/dth_it[jth],
      g0 = 1. - g1;

    const double
      h1 = (phi - phi_i)/dphi,
      h0 = 1. - h1;

    grid_data_point_t
      x000 = grid3d_data[IND3(it,jth,kphi)],
      x001 = grid3d_data[IND3(it,jth,kphi1)],
      x010 = grid3d_data[IND3(it,jth1,kphi)],
      x011 = grid3d_data[IND3(it,jth1,kphi1)],
      x100 = grid3d_data[IND3(it1,jth,kphi)],
      x101 = grid3d_data[IND3(it1,jth,kphi1)],
      x110 = grid3d_data[IND3(it1,jth1,kphi)],
      x111 = grid3d_data[IND3(it1,jth1,kphi1)];

    retval = f0*g0*h0*x000
           + f0*g0*h1*x001
           + f0*g1*h0*x010
           + f0*g1*h1*x011
           + f1*g0*h0*x100
           + f1*g0*h1*x101
           + f1*g1*h0*x110
           + f1*g1*h1*x111;

/*
if (retval.rho < 0) {
  using namespace std;
  log_one(error) << "internal: negative density!" << endl;
  cout << "negative density ("<< retval.rho << ")at:" << endl;
  cout << " - t     = " << tm    << ", t_i     = " << tm_i << endl;
  cout << " - theta = " << theta << ", theta_i = " << theta_i << endl;
  cout << " - phi   = " << phi   << ", phi_i   = " << phi_i << endl;
  cout << " - {it, jth, kphi} = " << it <<","<< jth<<","<< kphi << endl;
  cout << "x000: " << x000.rho << endl;
  cout << "x001: " << x001.rho << endl;
  cout << "x010: " << x010.rho << endl;
  cout << "x011: " << x011.rho << endl;
  cout << "x100: " << x100.rho << endl;
  cout << "x101: " << x101.rho << endl;
  cout << "x110: " << x110.rho << endl;
  cout << "x111: " << x111.rho << endl;
  exit (0);
}
*/

  }

  return retval;
}

/**
* @brief   Initial scan of the flux data and parameter setup
*
*/
void
init_read_ascii_flux_files() {
  using namespace std;
  using namespace param;
  string line;

  log_one(info) << "reading flux files at '"
                << input_flux_files << "'" << endl;
  glob_input_filenames(input_flux_files);
  sort(input_filenames.begin(), input_filenames.end());
  read_time_stamps();

  // consistency check: see that all extraction radii are identical
  extraction_radius = read_extraction_radius(input_filenames[0]);
  for(auto fname : input_filenames) {
    double R_ext = read_extraction_radius(fname);
    mpi_assert(R_ext == extraction_radius);
  }
  log_one(info) << "extraction radius: "
                << extraction_radius << " [cm]" << endl;

  // set grid size
  INFLX_NT = grid_times.size();
  auto [Nth, Nph] = read_spherical_grid_dimensions(input_filenames[0]);
  INFLX_NTHETA = Nth;
  INFLX_NPHI = Nph;

  // resize data arrays
  grid2d_theta.resize(INFLX_NT*INFLX_NTHETA);
  grid2d_dth.resize(INFLX_NT*INFLX_NTHETA);
  grid3d_data.resize(INFLX_NT*INFLX_NPHI*INFLX_NTHETA);
  grid3d_cumulative_mass.resize(INFLX_NT*INFLX_NPHI*INFLX_NTHETA);
  grid1d_cumulative_mass.resize(INFLX_NT);
  total_ejecta_mass = compute_total_mass();

  log_one(info) << "total mass of the injected flux: "
                << total_ejecta_mass << " [Msun]" << endl;

/*
// output cumulative mass and mass loss rate as a function of time
// TODO: output by request in a user-specified output file
for (int i = 1; i<INFLX_NT-1; ++i) {
  double t = 0.5*(grid_times[i] + grid_times[i-1]);
  double dt = grid_times[i] - grid_times[i-1];
  double m1 = grid1d_cumulative_mass[i-1];
  double m2 = grid1d_cumulative_mass[i];
  printf("%12.5e  %12.5e  %12.5e\n", t, m2, (m2-m1)/dt);
}
*/

//cout << "total mass: " << (mass/M_SUN_CGS) << endl;
//auto x = linear_interpolator(0.3, M_PI/2., M_PI);
//cout << "density: " << x.rho << endl;
//exit(0);
}

/**
* @brief   Reads the HDF5 flux file (which is actually a ballistically-expanded
*          density)
*
*/
void
init_read_hdf5_flux_file() {
    hid_t h5file;
    h5file = h5aux::H5P_openFile(param::input_flux_files, H5F_ACC_RDONLY);
    std::vector<std::string> datasets;
    
    // list the datasets in the HDF5 file
    herr_t status = H5Literate(h5file, H5_INDEX_NAME, H5_ITER_NATIVE, NULL, 
                    h5aux::list_datasets, &datasets);
    std::cout << "Datasets in '" << param::input_flux_files << "': ";
    for (const auto &name : datasets)
        std::cout << name << " ";
    std::cout << std::endl;

    int ndims;
    status = h5aux::H5D_getDimensions(h5file, "density", &ndims, grid3d_dims);
    std::cout << "Dimensions of the 'density' dataset: [" 
              << grid3d_dims[0] <<","
              << grid3d_dims[1] <<","
              << grid3d_dims[2] <<"]"<< std::endl;
  
    hsize_t grid3d_npts = grid3d_dims[0]*grid3d_dims[1]*grid3d_dims[2];
    double *data = new double[grid3d_npts];
    status = h5aux::H5D_readDataset(h5file, "density", data);
    std::cout << data[10345] << std::endl;
    H5Fclose(h5file);
    delete[] data;
    return; // TODO
}

/**
* @brief   Initialize the influx namespace
*
*/
void
init() {
  using namespace param;

  // If the last two characters of the `input_flux_files` is `*.h5`,
  // read a single HDF5 file; otherwise, read bunch of files
  size_t inpflen = strlen(input_flux_files);
  if (input_flux_files[inpflen-2] == 'h' && input_flux_files[inpflen-1] == '5') {
      log_one(info) << "reading HDF5 flux file at '"
                << input_flux_files << "'" << std::endl;
      init_read_hdf5_flux_file();
  }
  else {
      log_one(info) << "reading ASCII flux files at '"
                << input_flux_files << "'" << std::endl;
      init_read_ascii_flux_files();
  }

  //exit(0); // DEBUG
}

#undef IND3
#undef IND2

} // namespace influx


/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2017 Triad National Security, LLC
 * All rights reserved.
 *~--------------------------------------------------------------------------~*/

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
 * @file density_profiles.h
 * @author Oleg Korobkin
 * @date January 2019
 * @brief Interface to select various density profiles
 *
 * Notes:
 *  - all spherical density profiles have support with radius R = 1;
 *  - all profiles are normalized to total mass M = 1.
 */

#ifndef DENSITY_PROFILES_H
#define DENSITY_PROFILES_H

#include "params.h"
#include "tree.h"
#include "user.h"
#include "math.h"
#include <boost/algorithm/string.hpp>
#include <math.h>
#include <stdlib.h>
#include "lane_emden.h"
#include "h5aux.h"
#include <H5Cpp.h>

namespace density_profiles {

// spherical density profile function
typedef double (*radial_function_t)(const double);
static radial_function_t spherical_density_profile = NULL;
static radial_function_t spherical_mass_profile = NULL;
static radial_function_t spherical_drho_dr = NULL;
static radial_function_t spherical_alpha2 = NULL;
static radial_function_t spherical_dalpha2_dr = NULL;
static radial_function_t spherical_beta2 = NULL;
static radial_function_t spherical_dbeta2_dr = NULL;

// N-dimensional (N=2,3) density profile function and its gradient
typedef double (*function_ndim_t)(const point_t &);
typedef point_t (*grad_function_ndim_t)(const point_t &);
static function_ndim_t density_ndim = NULL;
static grad_function_ndim_t grad_density_ndim = NULL;

// constants for the mesa density
static double mesa_rho0;
static double mesa_q; // ratio of the slope width to the radius

// tabulated density profiles
static std::vector<double> rad_grid;
static std::vector<double> rho_grid;
static std::vector<double> mass_grid;
static std::vector<double> drhodr_grid;
static std::vector<double> alpha2_grid;
static std::vector<double> dalpha2dr_grid;
static std::vector<double> beta2_grid;
static std::vector<double> dbeta2dr_grid;

static std::vector<double> theta_grid;
static std::vector<double> phi_grid;
static double phi_grid_min;

static std::vector<double> p_grid;
static std::vector<double> ie_grid;
static std::vector<double> eps_grid;
static std::vector<double> ye_grid;
static std::vector<double> temp_grid;
static std::vector<double> vr_grid;
static std::vector<double> vt_grid;
static std::vector<double> vp_grid;

// components of the density gradient in 3D
static std::vector<double> drhodx_grid;
static std::vector<double> drhody_grid;
static std::vector<double> drhodz_grid;

// interpolator objects
//
static interp::linear_interpolator_3d_nug rho_interp;
static interp::linear_interpolator_3d_nug p_interp;
static interp::linear_interpolator_3d_nug ie_interp;
static interp::linear_interpolator_3d_nug eps_interp;
static interp::linear_interpolator_3d_nug temp_interp;
static interp::linear_interpolator_3d_nug ye_interp;
static interp::linear_interpolator_3d_nug vr_interp;
static interp::linear_interpolator_3d_nug vt_interp;
static interp::linear_interpolator_3d_nug vp_interp;
static interp::linear_interpolator_3d_nug drho_dx_interp;
static interp::linear_interpolator_3d_nug drho_dy_interp;
static interp::linear_interpolator_3d_nug drho_dz_interp;

double rho_ndim_from_data_grid(const point_t & rp);
point_t grad_rho_ndim_from_data_grid(const point_t & rp);

/**
 * @brief  constant uniform density in a domain of radius R = 1,
 *         normalized such that the total mass M = 1
 * @param  r     - spherical radius
 */
double
rho_constant_density(const double r) {
  double rho = 0.0;
  if constexpr(gdimension == 1)
    rho = 0.5;

  if constexpr(gdimension == 2)
    rho = 1.0 / M_PI;

  if constexpr(gdimension == 3)
    rho = 0.75 / M_PI;
  return rho;
}

double
mass_constant_density(const double r) {
  double mass = 0.0;
  if constexpr(gdimension == 1)
    mass = r;

  if constexpr(gdimension == 2)
    mass = SQ(r);

  if constexpr(gdimension == 3)
    mass = CU(r);

  return mass;
}

double
drhodr_constant_density(const double r) {
  return 0.;
}

/**
 * @brief  parabolic density
 * @param  r     - spherical radius
 */
double
rho_parabolic_density(const double r) {
  double rho = 0.0;
  if constexpr(gdimension == 1)
    rho = 0.75 * (1. - SQ(r));

  if constexpr(gdimension == 2)
    rho = 2. / M_PI * (1. - SQ(r));

  if constexpr(gdimension == 3)
    rho = 15. / (8. * M_PI) * (1. - SQ(r));

  return rho;
}

double
mass_parabolic_density(const double r) {
  double mass = 0.0;
  if constexpr(gdimension == 1)
    mass = 0.5 * r * (3.0 - SQ(r));

  if constexpr(gdimension == 2)
    mass = SQ(r) * (2. - SQ(r));

  if constexpr(gdimension == 3)
    mass = 0.5 * CU(r) * (5. - 3. * SQ(r));

  return mass;
}

double
drhodr_parabolic_density(const double r) {
  double drhodr = 0.0;
  if constexpr(gdimension == 1)
    drhodr = -1.5 * r;

  if constexpr(gdimension == 2)
    drhodr = -4. / M_PI * r;

  if constexpr(gdimension == 3)
    drhodr = -15. / (4. * M_PI) * r;

  return drhodr;
}

/**
 * @brief  spherical "mesa" density: flat top and steep slopes
 *
 *           / rho0                      if r < r0;
 *           |
 * rho(r) = <  rho0 (1 - (r-r0)^2/dr^2)  if r0 < r < 1;
 *           |
 *           \ 0                         if r > 1.
 *
 * @param  r     - spherical radius
 */
double
mesa_mass_helper(const double r) {
  const double dr = mesa_q, r0 = 1. - mesa_q;
  double mm = 0.0;

  if constexpr(gdimension == 2)
    mm = SQ(r) * (1. - (.5 * SQ(r) + SQ(r0)) / SQ(dr)) +
         (.5 * CU(r0) + 4. * CU(r)) * r0 / (3. * SQ(dr));

  if constexpr(gdimension == 3)
    mm = CU(r0) / 3. + (CU(r) - CU(r0)) / 3. * (1. - SQ(r0) / SQ(dr)) +
         r0 * (SQ(r) - SQ(r0)) * (SQ(r) + SQ(r0)) / (2. * SQ(dr)) -
         (SQ(r) * CU(r) - SQ(r0) * CU(r0)) / (5 * SQ(dr));
  return mm;
}

double
rho_mesa_density(const double r) {
  double rho = 0.0;
  const double r0 = 1. - mesa_q;
  if(r < 1. - mesa_q)
    rho = mesa_rho0;
  else if(r < .9999)
    rho = mesa_rho0 * (1. - SQ(r - r0) / SQ(mesa_q));
  return rho;
}

double
mass_mesa_density(const double r) {
  const double dr = mesa_q, r0 = 1. - mesa_q;
  double m = 0.0;
  if constexpr(gdimension == 1) {
    if(r < 1. - mesa_q)
      m = 2. * mesa_rho0 * r;
    else if(r - 1. < 1e-12)
      m = 2. * mesa_rho0 * (r - CU(r - r0) / (3. * SQ(dr)));
  }
  if constexpr(gdimension == 2) {
    if(r < 1. - mesa_q)
      m = M_PI * mesa_rho0 * SQ(r);
    else if(r - 1. < 1e-12)
      m = M_PI * mesa_rho0 * mesa_mass_helper(r);
  }
  if constexpr(gdimension == 3) {
    if(r < 1. - mesa_q)
      m = 4. * M_PI / 3. * mesa_rho0 * CU(r);
    else if(r - 1. < 1e-12)
      m = 4. * M_PI * mesa_rho0 * mesa_mass_helper(r);
  }
  return m;
}

double
drhodr_mesa_density(const double r) {
  double drhodr = 0.0;
  const double r0 = 1. - mesa_q;
  if(r > 1. - mesa_q and r < .9999)
    drhodr = -2. * mesa_rho0 * (r - r0) / SQ(mesa_q);
  return drhodr;
}

/**
 * @brief  kilonova spherical-ejecta density
 * @param  r     - spherical radius
 */
double
rho_kn_ejecta(const double r) {
  double z = 1. - r*r;
  double rho = z*z*z;
  if constexpr(gdimension == 1)
    rho *= 35./32.;

  if constexpr(gdimension == 2)
    rho *= 4. / M_PI;

  if constexpr(gdimension == 3)
    rho *= 315. / (64. * M_PI);

  return rho;
}

double
mass_kn_ejecta(const double r) {
  double r2 = r*r;
  double mass = 0.;
  if constexpr(gdimension == 1)
    mass = 35./16.*r*(1. - r2*(1. - r2*(.6 - r2/7.)));

  if constexpr(gdimension == 2)
    mass = 4.*r2*(1 - r2*(1.5 - r2*(1 - .25*r2)));

  if constexpr(gdimension == 3)
    mass = 315./16.*r*r2*(1./3. - r2*(.6 - r2*(3./7. - r2/9.)));

  return mass;
}

double
drhodr_kn_ejecta(const double r) {
  double z = 1. - r*r;
  double drhodr = -6.*r*z*z;
  if constexpr(gdimension == 1)
    drhodr *= 35./32.;

  if constexpr(gdimension == 2)
    drhodr *= 4. / M_PI;

  if constexpr(gdimension == 3)
    drhodr *= 315. / (64. * M_PI);

  return drhodr;
}

double
rho_ndim_kn_ejecta(const point_t & rp) {
  using namespace param;
  const double x = flecsi::magnitude(rp) / sphere_radius,
             rho = rho_initial*rho_kn_ejecta(x)/rho_kn_ejecta(0.);
  return rho;
}

point_t
grad_rho_kn_ejecta(const point_t & rp) {
  using namespace param;
  const double x = flecsi::magnitude(rp) / sphere_radius,
          drhodr = drhodr_kn_ejecta(x)*rho_initial/rho_kn_ejecta(0.)
                 / sphere_radius;
  point_t nr{0};
  nr = rp / (flecsi::magnitude(rp) + 1e-16);
  return drhodr * nr;
}

/**
 * @brief  Sharp density profile
 * @param  r     - spherical radius
 * @reference Simialar as Rosswog 1911.13093 Sec.3.1
 */
double
rho_sharp_spherical(const double r) {
  double rho = 0.0;
  double rho0 = 0.2;
  double drho = 0.8;
  if (r < 0.5){
    rho = rho0 + drho;
  } else {
    rho = rho0;
  }
  return rho;
}

double
mass_sharp_spherical(const double r) {
  double r2 = r*r;
  double mass = 0.;
  if constexpr(gdimension == 1)
    mass = r;

  if constexpr(gdimension == 2)
    mass = SQ(r)*M_PI;

  if constexpr(gdimension == 3)
    mass = 4./3.*CU(r)*M_PI;

  return mass;
}

double
drhodr_sharp_spherical(const double r) {
  return 0.;
}


/**
 * @brief  read the density input file
 * @param  ifname - 4-column ASCII file: 1:r 2:rho 3:m 4:drho/dr
 *                  possibly with a header with lines starting with '#'
 */
void
read_input_density_file(const char * ifname) {
  using namespace std;
  ifstream infile;
  string line;
  int ln, Nr;
  std::string::size_type sz1, sz2;
  double rad, rho, mass, drhodr;

  // attempt to open the file
  infile.open(ifname);
  if(!infile) {
    cerr << "ERROR: Unable to open density profile '" << ifname << "'" << endl;
    exit(1);
  }

  // count the number of lines, skipping the header
  Nr = 0;
  while(std::getline(infile, line)) {
    if(line.find("#") != string::npos or
       line.find_first_not_of(' ') == string::npos)
      continue;
    ++Nr;
  }
  cout << "Read density profile " << ifname << " with " << Nr << " data points."
       << endl;

  // allocate arrays
  rad_grid.resize(Nr);
  rho_grid.resize(Nr);
  mass_grid.resize(Nr);
  drhodr_grid.resize(Nr);

  // read in the values
  infile.clear();
  infile.seekg(0, ios::beg);
  for(ln = 0; std::getline(infile, line);) {
    if(line.find("#") != string::npos or
       line.find_first_not_of(' ') == string::npos)
      continue;
    rad_grid[ln] = std::stod(line, &sz1);
    rho_grid[ln] = std::stof(line.substr(sz1), &sz2);
    sz1 += sz2;
    mass_grid[ln] = std::stof(line.substr(sz1), &sz2);
    sz1 += sz2;
    drhodr_grid[ln] = std::stof(line.substr(sz1));

    ln++;
  }
}


/**
 * @brief  read the density input file
 * @param  ifname - hdf5 file: 3 dimensional array with density field
 *             and spherical grid coordinates
 */
void
read_input_density_h5file(const char * ifname) {
   static bool called_once = false;
   if (called_once) return;
   called_once = true;

   //step 0: open file
   hid_t file_id = h5aux::H5P_openFile( ifname, H5F_ACC_RDONLY);
   //step 1: read meta data in file, group density
   hsize_t dims[3];
   int ndims;
   h5aux::H5D_getDimensions(file_id, "density", &ndims, dims);
   //step 2: read coordinate grid data
   rad_grid.resize(dims[2]);
   h5aux::H5D_readDataset(file_id, "r",&(rad_grid[0]));
   
   theta_grid.resize(dims[0] + 2);
   h5aux::H5D_readDataset(file_id, "theta",&(theta_grid[1]));

   phi_grid.resize(dims[1] + 1);
   h5aux::H5D_readDataset(file_id, "phi",&(phi_grid[0]));
   phi_grid_min = phi_grid[0];
   //step 3: allocate memory for density
   std::vector<double> rho_tmp;
   std::vector<double> p_tmp;
   std::vector<double> ie_tmp;
   std::vector<double> eps_tmp;
   std::vector<double> temp_tmp;
   std::vector<double> ye_tmp;
   std::vector<double> vr_tmp;
   std::vector<double> vt_tmp;
   std::vector<double> vp_tmp;
   rho_tmp.resize(dims[0]*dims[1]*dims[2]);
   h5aux::H5D_read3DDataset(file_id, "density",&(rho_tmp[0]));

   p_tmp.resize(dims[0]*dims[1]*dims[2]);
   h5aux::H5D_read3DDataset(file_id, "pressure",&(p_tmp[0]));
   
   ie_tmp.resize(dims[0]*dims[1]*dims[2]);
   h5aux::H5D_read3DDataset(file_id, "int_energy",&(ie_tmp[0]));

   eps_tmp.resize(dims[0]*dims[1]*dims[2]);
   h5aux::H5D_read3DDataset(file_id, "eps",&(eps_tmp[0]));
   
   temp_tmp.resize(dims[0]*dims[1]*dims[2]);
   h5aux::H5D_read3DDataset(file_id, "temperature",&(temp_tmp[0]));

   ye_tmp.resize(dims[0]*dims[1]*dims[2]);
   h5aux::H5D_read3DDataset(file_id, "y_e",&(ye_tmp[0]));

   vr_tmp.resize(dims[0]*dims[1]*dims[2]);
   h5aux::H5D_read3DDataset(file_id, "radial_vel",&(vr_tmp[0]));

   vt_tmp.resize(dims[0]*dims[1]*dims[2]);
   h5aux::H5D_read3DDataset(file_id, "theta_vel",&(vt_tmp[0]));

   vp_tmp.resize(dims[0]*dims[1]*dims[2]);
   h5aux::H5D_read3DDataset(file_id, "phi_vel",&(vp_tmp[0]));
   
   //step 4: read density field
   rho_grid.resize((dims[0] + 2)*(dims[1] + 1)*dims[2]);
   p_grid.resize((dims[0] + 2)*(dims[1] + 1)*dims[2]);
   ie_grid.resize((dims[0] + 2)*(dims[1] + 1)*dims[2]);
   eps_grid.resize((dims[0] + 2)*(dims[1] + 1)*dims[2]);
   temp_grid.resize((dims[0] + 2)*(dims[1] + 1)*dims[2]);
   ye_grid.resize((dims[0] + 2)*(dims[1] + 1)*dims[2]);
   vr_grid.resize((dims[0] + 2)*(dims[1] + 1)*dims[2]);
   vt_grid.resize((dims[0] + 2)*(dims[1] + 1)*dims[2]);
   vp_grid.resize((dims[0] + 2)*(dims[1] + 1)*dims[2]);
   for (int kth = 0; kth < dims[0]; kth++) {
       for (int ir = 0; ir < dims[2]; ir++) {
           for (int jph = 0; jph < dims[1]; jph++) {
               int ijk  = ir + dims[2]*(jph + dims[1]*kth);
               int ijk1 = ir + dims[2]*(jph + (dims[1] + 1)*(kth + 1));
               rho_grid[ijk1] = rho_tmp[ijk];
               p_grid[ijk1] = p_tmp[ijk];
               ie_grid[ijk1] = ie_tmp[ijk];
               eps_grid[ijk1] = eps_tmp[ijk];
               temp_grid[ijk1] = temp_tmp[ijk];
               ye_grid[ijk1] = ye_tmp[ijk];
               vr_grid[ijk1] = vr_tmp[ijk];
               vt_grid[ijk1] = vt_tmp[ijk];
               vp_grid[ijk1] = vp_tmp[ijk];
           }
           int ijk1 = ir + dims[2]*(dims[1] + (dims[1] + 1)*(kth + 1));
           int ijk0 = ir + dims[2]*(0       + (dims[1] + 1)*(kth + 1));
           rho_grid[ijk1] = rho_grid[ijk0];
           p_grid[ijk1] = p_grid[ijk0];
           ie_grid[ijk1] = ie_grid[ijk0];
           eps_grid[ijk1] = eps_grid[ijk0];
           temp_grid[ijk1] = temp_grid[ijk0];
           ye_grid[ijk1] = ye_grid[ijk0];
           vr_grid[ijk1] = vr_grid[ijk0];
           vt_grid[ijk1] = vt_grid[ijk0];
           vp_grid[ijk1] = vp_grid[ijk0];
       }
   }

   for (int ir = 0; ir < dims[2]; ir++) {
       for (int jph = 0; jph < dims[1] + 1; jph++) {
           int jpha = (jph + dims[1]/2) % dims[1];
           int ijk1 = ir + dims[2]*(jpha + (dims[1] + 1)*1);
           int ijk0 = ir + dims[2]*(jph  + (dims[1] + 1)*0);
           rho_grid[ijk0] = rho_grid[ijk1];

           ijk1 = ir + dims[2]*(jpha + (dims[1] + 1)*(dims[0]));
           ijk0 = ir + dims[2]*(jph  + (dims[1] + 1)*(dims[0] + 1));
           rho_grid[ijk0] = rho_grid[ijk1];
           p_grid[ijk0] = p_grid[ijk1];
           ie_grid[ijk0] = ie_grid[ijk1];
           eps_grid[ijk0] = eps_grid[ijk1];
           temp_grid[ijk0] = temp_grid[ijk1];
           ye_grid[ijk0] = ye_grid[ijk1];
           vr_grid[ijk0] = vr_grid[ijk1];
           vt_grid[ijk0] = vt_grid[ijk1];
           vp_grid[ijk0] = vp_grid[ijk1];
       }
   }
   phi_grid[dims[1]++] = phi_grid[0] + 2*M_PI;
   theta_grid[0] = -theta_grid[1];
   theta_grid[dims[0] + 1] = 2*M_PI - theta_grid[dims[0]];
   dims[0] += 2;

  // printf("the first  element of the theta grid is %e \n", theta_grid[0] * 180 / M_PI);
  // printf("the second element of the theta grid is %e \n", theta_grid[1] * 180 / M_PI);
  // printf("the last element of the theta grid is %e \n", theta_grid[dims[0] - 1] * 180 / M_PI);

   // initialize our interpolators
   rho_interp.set_data(&(rad_grid[0]),&(phi_grid[0]),&(theta_grid[0]),
           &(rho_grid[0]),dims[2],dims[1],dims[0]);

   p_interp.set_data(&(rad_grid[0]),&(phi_grid[0]),&(theta_grid[0]),
           &(p_grid[0]),dims[2],dims[1],dims[0]);
   ie_interp.set_data(&(rad_grid[0]),&(phi_grid[0]),&(theta_grid[0]),
           &(ie_grid[0]),dims[2],dims[1],dims[0]);
   eps_interp.set_data(&(rad_grid[0]),&(phi_grid[0]),&(theta_grid[0]),
           &(eps_grid[0]),dims[2],dims[1],dims[0]);
   temp_interp.set_data(&(rad_grid[0]),&(phi_grid[0]),&(theta_grid[0]),
           &(temp_grid[0]),dims[2],dims[1],dims[0]);
   ye_interp.set_data(&(rad_grid[0]),&(phi_grid[0]),&(theta_grid[0]),
           &(ye_grid[0]),dims[2],dims[1],dims[0]);
   vr_interp.set_data(&(rad_grid[0]),&(phi_grid[0]),&(theta_grid[0]),
           &(vr_grid[0]),dims[2],dims[1],dims[0]);
   vt_interp.set_data(&(rad_grid[0]),&(phi_grid[0]),&(theta_grid[0]),
           &(vt_grid[0]),dims[2],dims[1],dims[0]);
   vp_interp.set_data(&(rad_grid[0]),&(phi_grid[0]),&(theta_grid[0]),
           &(vp_grid[0]),dims[2],dims[1],dims[0]);

//for (double rr=1e9; rr<1e10; rr+=1e9) {
//printf("vr_interp(%e,phi_grid[15],theta_grid[0]..[1].....[-2]..[-1])= %e  %e     %e  %e\n",
//    rr,
//    vr_interp(rr,phi_grid[15],theta_grid[0]),
//    vr_interp(rr,phi_grid[15],theta_grid[1]),
//    vr_interp(rr,phi_grid[15],theta_grid[dims[0]-2]),
//    vr_interp(rr,phi_grid[15],theta_grid[dims[0]-1])
//    );
//
//}
//for (double th=0; th<M_PI+0.1; th+=M_PI/20) {
//printf("vr_interp(rad_grid[400],phi_grid[1/65], %fpi)  %e  %e\n",
//    th/M_PI,
//    vr_interp(rad_grid[400],phi_grid[1],th),
//    vr_interp(rad_grid[400],phi_grid[dims[1]-1],th)
//    );
//}
   //point_t rp {2e9,2e9,3e9};
   //printf("rho_ndim_from_data_grid()=%e \n", rho_ndim_from_data_grid(rp));

   //step 6: compute density gradients in spherical coords
   //        we need to know the index order
   //        example: datasets have shape {128, 66, 999}
   //                 - index for theta runs from 0 .. 127
   //                 - index for phi   runs from 0 .. 65
   //                 - index for r     runs from 0 .. 998;
   //        index for r is the fastest-changing index
   //

   std::vector<double> drhodtheta_grid;
   std::vector<double> drhodphi_grid;

   drhodr_grid.resize(dims[0]*dims[1]*dims[2]);
   drhodtheta_grid.resize(dims[0]*dims[1]*dims[2]);
   drhodphi_grid.resize(dims[0]*dims[1]*dims[2]);

   std::vector<double> xp, yp, zp;
   xp.resize(dims[0]*dims[1]*dims[2]);
   yp.resize(dims[0]*dims[1]*dims[2]);
   zp.resize(dims[0]*dims[1]*dims[2]);

   for (int kth = 0; kth < dims[0]; kth++) {
       int kp = std::min(kth + 1, (int)(dims[0]-1));
       int km = std::max(kth - 1, 0);
       double dtheta = theta_grid[kp] - theta_grid[km];
       for (int jph = 0; jph < dims[1]; jph++) {
           int jp = jph + 1;
           int jm = jph ? (jph - 1) : (dims[1] - 2);
           double dphi = phi_grid[jp] - phi_grid[jm];
           if (jph == 0) dphi = phi_grid[dims[1]-1] - phi_grid[dims[1]-2];
           for (int ir = 0; ir < dims[2]; ir++) {
               int ijk = ir + dims[2]*(jph + dims[1]*kth);
               int ip  = std::min(ir + 1, (int)(dims[2]-1));
               int im  = std::max(ir - 1, 0);
               double dr = rad_grid[ip] - rad_grid[im];

               int ijkp = ip + dims[2]*(jph + dims[1]*kth);
               int ijkm = im + dims[2]*(jph + dims[1]*kth);
               drhodr_grid[ijk] = (rho_grid[ijkp] - rho_grid[ijkm])/dr;

               ijkp = ir + dims[2]*(jph + dims[1]*kp);
               ijkm = ir + dims[2]*(jph + dims[1]*km);
               drhodtheta_grid[ijk] = (rho_grid[ijkp] - rho_grid[ijkm])/dtheta;

               ijkp = ir + dims[2]*(jp + dims[1]*kth);
               ijkm = ir + dims[2]*(jm + dims[1]*kth);
               drhodphi_grid[ijk] = (rho_grid[ijkp] - rho_grid[ijkm])/dphi;
//if (jph == 30 && kth == 20) {
//printf("%d  %12.5f  %14.7e %14.7e\n", ir, rad_grid[ir], rho_grid[ijk],  drhodr_grid[ijk]);
//}
           }
       }
   }

   for (int ir = 0; ir < dims[2]; ir++) {
       for (int jph = 0; jph < dims[1]; jph++) {
           int jpha = (jph + dims[1]/2) % dims[1];
           int ijk1 = ir + dims[2]*(jpha + dims[1]*1);
           int ijk0 = ir + dims[2]*(jph  + dims[1]*0);
           drhodr_grid[ijk0]     = drhodr_grid[ijk1];
           drhodtheta_grid[ijk0] = drhodtheta_grid[ijk1];
           drhodphi_grid[ijk0]   = drhodphi_grid[ijk1];

           ijk1 = ir + dims[2]*(jpha + dims[1]*(dims[0]-2));
           ijk0 = ir + dims[2]*(jph  + dims[1]*(dims[0]-1));
           drhodr_grid[ijk0]     = drhodr_grid[ijk1];
           drhodtheta_grid[ijk0] = drhodtheta_grid[ijk1];
           drhodphi_grid[ijk0]   = drhodphi_grid[ijk1];
       }
   }


   //step 7: transform these gradients to cartesian frame
   drhodx_grid.resize(dims[0]*dims[1]*dims[2]);
   drhody_grid.resize(dims[0]*dims[1]*dims[2]);
   drhodz_grid.resize(dims[0]*dims[1]*dims[2]);
   for (int kth = 0; kth < dims[0]; kth++) {
       double th = theta_grid[kth],
           costh = cos(th),
           sinth = sin(th);
       for (int jph = 0; jph < dims[1]; jph++) {
           double phi = phi_grid[jph],
               cosphi = cos(phi),
               sinphi = sin(phi);
           for (int ir = 0; ir < dims[2]; ir++) {
               int ijk = ir + dims[2]*(jph + dims[1]*kth);
               double eps = 1e-15*param::sphere_radius;
               double r = rad_grid[ir];
               double z = r*costh;
               double R = r*sinth;
               double x = R*cosphi;
               double y = R*sinphi;

               xp[ijk] = x;
               yp[ijk] = y;
               zp[ijk] = z;

               // do the transformation!
               double drho_dr_pt = drhodr_grid[ijk];
               double drho_dph_pt = drhodphi_grid[ijk];
               double drho_dth_pt = drhodtheta_grid[ijk];

               // drhodx = drdx*drhodr + dthdx*drhodth + dphidx*drhodphi
               drhodx_grid[ijk] = x/(r + eps)       *drho_dr_pt
                                + x*z/(r*r*R + eps) *drho_dth_pt
                                - y/(R*R + eps)     *drho_dph_pt;

               // drhody = drdy*drhodr + dthdy*drhodth + dphidy*drhodphi
               drhody_grid[ijk] = y/(r + eps)       *drho_dr_pt
                                + y*z/(r*r*R + eps) *drho_dth_pt
                                + x/(R*R + eps)     *drho_dph_pt;

               // drhodz = drdz*drhodr + dthdz*drhodth + dphidz*drhodphi
               drhodz_grid[ijk] = z/(r + eps)       *drho_dr_pt
                                - R/(r*r + eps)     *drho_dth_pt;

//if (jph == 30 && kth == 20) {
//printf("%d  %12.5f  %14.7e %14.7e\n", ir, rad_grid[ir], rho_grid[ijk],  drhodz_grid[ijk]);
//}
           }
       }
   }

   // setup the interpolators for grad_rho
   drho_dx_interp.set_data(&(rad_grid[0]),&(phi_grid[0]),&(theta_grid[0]),
           &(drhodx_grid[0]),dims[2],dims[1],dims[0]);
   drho_dy_interp.set_data(&(rad_grid[0]),&(phi_grid[0]),&(theta_grid[0]),
           &(drhody_grid[0]),dims[2],dims[1],dims[0]);
   drho_dz_interp.set_data(&(rad_grid[0]),&(phi_grid[0]),&(theta_grid[0]),
           &(drhodz_grid[0]),dims[2],dims[1],dims[0]);

   //step 8: close group and close file
   H5Fclose(file_id);

   //step 9: create file to test gradient computations
   // Reusing code from app/id_generators/regrid/main.cc
   file_id = H5Fcreate("gradient_test.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
   hid_t     dataset_id, dataspace_id, status;
   // create dataspace / dataset

   // record the coordinates
   dataspace_id = H5Screate_simple(3, dims, NULL);

   // x
   dataset_id = H5Dcreate(file_id, "x", H5T_NATIVE_DOUBLE,
                          dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
   status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,
                     H5S_ALL, H5P_DEFAULT, xp.data());
   status = H5Dclose(dataset_id);

   // y
   dataset_id = H5Dcreate(file_id, "y", H5T_NATIVE_DOUBLE,
                          dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
   status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,
                     H5S_ALL, H5P_DEFAULT, yp.data());
   status = H5Dclose(dataset_id);

   // z
   dataset_id = H5Dcreate(file_id, "z", H5T_NATIVE_DOUBLE,
                          dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
   status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,
                     H5S_ALL, H5P_DEFAULT, zp.data());
   status = H5Dclose(dataset_id);

   // record the density
   dataset_id = H5Dcreate(file_id, "rho", H5T_NATIVE_DOUBLE,
                          dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
   status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,
                     H5S_ALL, H5P_DEFAULT, rho_grid.data());
   status = H5Dclose(dataset_id);

   // record the gradients
   dataset_id = H5Dcreate(file_id, "grad_x", H5T_NATIVE_DOUBLE,
                          dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
   status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,
                     H5S_ALL, H5P_DEFAULT, drhodx_grid.data());
   status = H5Dclose(dataset_id);

   dataset_id = H5Dcreate(file_id, "grad_y", H5T_NATIVE_DOUBLE,
                          dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
   status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,
                     H5S_ALL, H5P_DEFAULT, drhody_grid.data());
   status = H5Dclose(dataset_id);

   dataset_id = H5Dcreate(file_id, "grad_z", H5T_NATIVE_DOUBLE,
                          dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
   status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,
                     H5S_ALL, H5P_DEFAULT, drhodz_grid.data());
   status = H5Dclose(dataset_id);

   status = H5Sclose(dataspace_id);
   status = H5Fclose(file_id);

   // Test gradient interpolation
   dims[0] = 140;
   dims[1] = 140;
   dims[2] = 140;
   double cube_side = param::sphere_radius;
   double dx = 2*cube_side / dims[0];
   static std::vector<double> rho_gr, drhodx_gr, drhody_gr, drhodz_gr;

   xp.resize(dims[0]*dims[1]*dims[2]);
   yp.resize(dims[0]*dims[1]*dims[2]);
   zp.resize(dims[0]*dims[1]*dims[2]);

   rho_gr.resize(dims[0]*dims[1]*dims[2]);
   drhodx_gr.resize(dims[0]*dims[1]*dims[2]);
   drhody_gr.resize(dims[0]*dims[1]*dims[2]);
   drhodz_gr.resize(dims[0]*dims[1]*dims[2]);

   for (int k = 0; k < dims[0]; k++) {
       double z = -cube_side + dx*k;
       for (int j = 0; j < dims[1]; j++) {
           double y = -cube_side + dx*j;
           for (int i = 0; i < dims[2]; i++) {
               double x = -cube_side + dx*i;
               int ijk = i + dims[2]*(j + dims[1]*k);

               xp[ijk] = x;
               yp[ijk] = y;
               zp[ijk] = z;

               point_t rp{0}, grad_rho{0};
               rp[0] = x;
               rp[1] = y;
               rp[2] = z;
               grad_rho = grad_rho_ndim_from_data_grid(rp);

               rho_gr[ijk]    = rho_ndim_from_data_grid(rp);
               drhodx_gr[ijk] = grad_rho[0];
               drhody_gr[ijk] = grad_rho[1];
               drhodz_gr[ijk] = grad_rho[2];
           }
       }
   }

   // Reusing code from app/id_generators/regrid/main.cc
   file_id = H5Fcreate("interp_gradient_test.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
   // create dataspace / dataset

   // record the coordinates
   dataspace_id = H5Screate_simple(3, dims, NULL);

   // x
   dataset_id = H5Dcreate(file_id, "x", H5T_NATIVE_DOUBLE,
                          dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
   status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,
                     H5S_ALL, H5P_DEFAULT, xp.data());
   status = H5Dclose(dataset_id);

   // y
   dataset_id = H5Dcreate(file_id, "y", H5T_NATIVE_DOUBLE,
                          dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
   status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,
                     H5S_ALL, H5P_DEFAULT, yp.data());
   status = H5Dclose(dataset_id);

   // z
   dataset_id = H5Dcreate(file_id, "z", H5T_NATIVE_DOUBLE,
                          dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
   status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,
                     H5S_ALL, H5P_DEFAULT, zp.data());
   status = H5Dclose(dataset_id);

   // record the density
   dataset_id = H5Dcreate(file_id, "rho", H5T_NATIVE_DOUBLE,
                          dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
   status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,
                     H5S_ALL, H5P_DEFAULT, rho_gr.data());
   status = H5Dclose(dataset_id);

   // record the gradients
   dataset_id = H5Dcreate(file_id, "grad_x", H5T_NATIVE_DOUBLE,
                          dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
   status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,
                     H5S_ALL, H5P_DEFAULT, drhodx_gr.data());
   status = H5Dclose(dataset_id);

   dataset_id = H5Dcreate(file_id, "grad_y", H5T_NATIVE_DOUBLE,
                          dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
   status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,
                     H5S_ALL, H5P_DEFAULT, drhody_gr.data());
   status = H5Dclose(dataset_id);

   dataset_id = H5Dcreate(file_id, "grad_z", H5T_NATIVE_DOUBLE,
                          dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
   status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,
                     H5S_ALL, H5P_DEFAULT, drhodz_gr.data());
   status = H5Dclose(dataset_id);

   status = H5Sclose(dataspace_id);
   status = H5Fclose(file_id);
}


/**
 * @brief   get index i such that xp[i] < x < xp[i+1] (binary search)
 * @param   x     - the value to localize;
 * @param   xp    - array of increasing values where to localize x.
 * @return  index i of an interval containing point x;
 *          '-1'   if x < x[0];
 *          'N-1'  if x > x[N-1] (N is the vector size).
 */
int
get_interval_index(const double x, const std::vector<double> & xp) {
  const int N = xp.size();
  if(x * (1 + 1e-15) < xp[0])
    return -1;
  if(x * (1 - 1e-15) > xp[N - 1])
    return N - 1;
  int i, i1 = 0, i2 = N - 1;
  while(i2 - i1 > 1) {
    i = (i1 + i2) / 2;
    if(x < xp[i])
      i2 = i;
    else
      i1 = i;
  }
  return i;
}

/**
 * @brief   cubic interpolation
 * @param   x        - the x-coordinate location where to interpolate
 * @param   xp       - the x-coordinates of data points: increasing
 * @param   yp       - the y-coordinates of data points
 * @return  if x is inside the range of xp, returns interpolated value;
 *          if x is outside the range of xp: returns zero
 */
double
cubic_interp(const double x,
  const std::vector<double> & xp,
  const std::vector<double> & yp) {
  int i2 = get_interval_index(x, xp);
  const int N = xp.size();
  if(i2 == 0)
    i2 = 1;
  if(i2 > N - 3)
    i2 = N - 3;
  int i1 = i2 - 1, i3 = i2 + 1, i4 = i2 + 2;

  double xx1 = x - xp[i1], xx2 = x - xp[i2], xx3 = x - xp[i3], xx4 = x - xp[i4];

  return yp[i1] * xx2 / (xx2 - xx1) * xx3 / (xx3 - xx1) * xx4 / (xx4 - xx1) +
         yp[i2] * xx1 / (xx1 - xx2) * xx3 / (xx3 - xx2) * xx4 / (xx4 - xx2) +
         yp[i3] * xx1 / (xx1 - xx3) * xx2 / (xx2 - xx3) * xx4 / (xx4 - xx3) +
         yp[i4] * xx1 / (xx1 - xx4) * xx2 / (xx2 - xx4) * xx3 / (xx3 - xx4);
}

/**
 * @brief  density profile from file, specified by
 *         the parameter density_profile_input
 * @param  r     - spherical radius
 */
double
rho_from_data_grid(const double r) {
  return cubic_interp(r, rad_grid, rho_grid);
}

double
mass_from_data_grid(const double r) {
  return cubic_interp(r, rad_grid, mass_grid);
}

double
drhodr_from_data_grid(const double r) {
  return cubic_interp(r, rad_grid, drhodr_grid);
}

double
alpha2_from_data_grid(const double r) {
  return cubic_interp(r,rad_grid, alpha2_grid);
}

double
dalpha2dr_from_data_grid(const double r) {
  return cubic_interp(r,rad_grid, dalpha2dr_grid);
}

double
beta2_from_data_grid(const double r) {
  return cubic_interp(r,rad_grid, beta2_grid);
}

double
dbeta2dr_from_data_grid(const double r) {
  return cubic_interp(r,rad_grid, dbeta2dr_grid);
}

double
rho_ndim_from_data_grid(const point_t & rp) {
  //
  double x=rp[0], y=rp[1], z=rp[2];
  double r = sqrt(x*x + y*y + z*z);
  double theta = atan2(sqrt(x*x+y*y),z);
  double phi = atan2(y,x);
  if (phi < phi_grid_min) phi += 2*M_PI;

  double rho = rho_interp(r, phi, theta);

  return rho;
}

double
Q_ndim_from_data_grid(const point_t & rp, interp::linear_interpolator_3d_nug& Q_interp) {
  //
  double x=rp[0], y=rp[1], z=rp[2];
  double r = sqrt(x*x + y*y + z*z);
  double theta = atan2(sqrt(x*x+y*y),z);
  double phi = atan2(y,x);
  if (phi < phi_grid_min) phi += 2*M_PI;

  return Q_interp(r, phi, theta);
}

point_t
grad_rho_ndim_from_data_grid(const point_t & rp) {
  //
  point_t grad{0};

  double x=rp[0], y=rp[1], z=rp[2];
  double r = sqrt(x*x + y*y + z*z);
  double theta = atan2(sqrt(x*x+y*y),z);
  double phi = atan2(y,x);
  if (phi < phi_grid_min) phi += 2*M_PI;

  grad[0] = drho_dx_interp(r,phi,theta);
  grad[1] = drho_dy_interp(r,phi,theta);
  grad[2] = drho_dz_interp(r,phi,theta);

  return grad;
}

/**
 * @brief      Density profile selector
 */
void
select() {
  using namespace param;
  std::string str_profile{density_profile};
  for(int c = 0; c < str_profile.length(); ++c)
    if(str_profile[c] == ' ')
      str_profile[c] = '_';
    else if(str_profile[c] == '-')
      str_profile[c] = '_';

  if(boost::iequals(str_profile, "constant")) {
    spherical_density_profile = rho_constant_density;
    spherical_mass_profile = mass_constant_density;
    spherical_drho_dr = drhodr_constant_density;
  }
  else if(boost::iequals(str_profile, "parabolic")) {
    spherical_density_profile = rho_parabolic_density;
    spherical_mass_profile = mass_parabolic_density;
    spherical_drho_dr = drhodr_parabolic_density;
  }
  else if(boost::iequals(str_profile, "mesa")) {
    spherical_density_profile = rho_mesa_density;
    spherical_mass_profile = mass_mesa_density;
    spherical_drho_dr = drhodr_mesa_density;
    mesa_q = mesa_rim_width;
    if constexpr(gdimension == 1)
      mesa_rho0 = .5 / (1. - mesa_q / 3.);

    if constexpr(gdimension == 2)
      mesa_rho0 = 1. / (M_PI * mesa_mass_helper(1.));

    if constexpr(gdimension == 3)
      mesa_rho0 = 1. / (4. * M_PI * mesa_mass_helper(1.));
  }
  else if(boost::iequals(str_profile, "kn_ejecta")) {
    spherical_density_profile = rho_kn_ejecta;
    spherical_mass_profile = mass_kn_ejecta;
    spherical_drho_dr = drhodr_kn_ejecta;
/*
for (double x = 0; x < 1.0; x += 0.01) {
  printf ("%6.2f  %13.7e  %13.7e  %14.7e\n", x,
      spherical_density_profile(x),
      spherical_mass_profile(x),
      spherical_drho_dr(x));
}
exit(0);
*/
  }
  else if(boost::iequals(str_profile, "kn_ejecta_3d")) {
    density_ndim = rho_ndim_kn_ejecta;
    grad_density_ndim = grad_rho_kn_ejecta;
  }
  else if(boost::iequals(str_profile, "sharp_spherical")) {
    spherical_density_profile = rho_sharp_spherical;
    spherical_mass_profile = mass_sharp_spherical;
    spherical_drho_dr = drhodr_sharp_spherical;
  }
  else if(boost::iequals(str_profile, "from_file")) {
    // read rho input file
    // - if the file is a text file (*.dat), assume it's 1D;
    // - if it's *.h5, assume it's a 3D density
    int l = strlen(input_density_file);
    if(input_density_file[l-2]=='h' && input_density_file[l-1]=='5') {
        read_input_density_h5file(input_density_file);
        density_ndim = rho_ndim_from_data_grid;
        grad_density_ndim = grad_rho_ndim_from_data_grid;
    }
    else {
        read_input_density_file(input_density_file);
        spherical_density_profile = rho_from_data_grid;
        spherical_mass_profile = mass_from_data_grid;
        spherical_drho_dr = drhodr_from_data_grid;
    }
//for (double x = 0; x < 1.0; x += 0.01) {
//  printf ("%6.2f  %13.7e  %13.7e  %14.7e\n", x,
//      density_ndim(x),
//      grad_density_ndim(x));
//}
//exit(0);
  }
  else if(boost::iequals(str_profile, "lane_emden")) {
    int N_r = lane_emden_radial_N;
    // invoke Lane-Emden solver to compute the profile on the fly
    lane_emden::solve(N_r, rad_grid, rho_grid, mass_grid, drhodr_grid,
                      alpha2_grid, dalpha2dr_grid, beta2_grid, dbeta2dr_grid);
    spherical_density_profile = rho_from_data_grid;
    spherical_mass_profile = mass_from_data_grid;
    spherical_drho_dr = drhodr_from_data_grid;
    spherical_alpha2 = alpha2_from_data_grid;
    spherical_dalpha2_dr = dalpha2dr_from_data_grid;
    spherical_beta2 = beta2_from_data_grid;
    spherical_dbeta2_dr = dbeta2dr_from_data_grid;
  }
  else {
    logm(error) << "ERROR: wrong parameter in density_profiles";
    exit(2);
  }

} // select()

} // namespace density_profiles

#endif // DENSITY_PROFILES_H

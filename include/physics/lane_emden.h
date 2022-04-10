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
 * @file lane_emden.h
 * @brief Lane-Emden solver for spherical stellar structure
 *
 * This file contains a Lane-Emden solver for predefined equation of state.
 * Instead of integrating density (rho) and mass (m) along the radius r, it
 * integrates variables s:= r^2 and mass along the variable theta, where
 * theta is defined normalized density: rho = rho_c * theta^n. Such inverse
 * integration allows both ends of the integrated variable to be defined:
 * theta \in [1.0, 0.0]. The use of s = r**2 allows to exploit an analytic
 * approximate solution for the initial step, where r = 0.
 *
 * The code shows convergence rate = 1.5 (TODO: check again)
 *
 * @author Bing-Jyun Tsao, Oleg Korobkin
 * @date September 2020
 */
#pragma once

#include <iostream>
#include <cmath>
#include <vector>
#include <boost/algorithm/string.hpp>
#include <fstream>
#include <cstdio>
#include "body.h"
#include "eos.h"
#include "params.h"
#include "units.h"
namespace lane_emden {

/**
* @brief      Returns derivatives dm/dtheta and ds/dtheta
*
* @param      m      enclosed mass
* @param      s      s := r^2: squared radius
* @param      th     normalized density coordinate (rho = rho_c theta^n)
* @param      rho_c  central density
* @param      n      polytropic index at the center of the star
* @param      pt     particle carrying the rest passive thermodynamic props
*/
std::pair<double, double>
dms_dth(const double m, const double s, const double th,
    const double rho_c, const double n, body & pt) {
  using namespace eos;
  double rho = rho_c * pow(th, n);
  pt.setDensity(rho);
  compute_pressure(pt);
  compute_soundspeed(pt);
  compute_internal_energy(pt);
  double p = pt.getPressure();
  double u = pt.getInternalenergy();
  double cs = pt.getSoundspeed();
  const double CLIGHT2 = phys::clight * phys::clight;
  double dPdrho = param::lane_emden_isothermal
                ? get_dpdrho_at_temp(pt)
                : cs*cs;
  // tov correction terms
  double GR_cor_ds = 1.0;
  double GR_cor_dm = 1.0;
  if(param::tov_correction){
    double GR_cor_ds1 = (1 + (rho*u + p) / (CLIGHT2*rho));
    double GR_cor_ds2 = (1 + (4*M_PI*sqrt(s*s*s)*p)/(m*CLIGHT2));
    double GR_cor_ds3 = (1 - (2*phys::GN*m)/(sqrt(s)*CLIGHT2));

    GR_cor_ds = GR_cor_ds3 / (GR_cor_ds1 * GR_cor_ds2);
    GR_cor_dm = (1 + u/CLIGHT2);
  }

  double dsdth = -2*n*sqrt(s*s*s)/(phys::GN * m * th) * dPdrho * GR_cor_ds;
  double dmdth = dsdth * 2*M_PI*sqrt(s)*rho * GR_cor_dm;
  return {dmdth, dsdth};
}

/**
* @brief      Runge-Kutta 4th-order integrator: single step
*
* @param      m      enclosed mass
* @param      s      s := r^2: squared radius
* @param      th     normalized density coordinate (rho = rho_c theta^n)
* @param      dth    theta step
* @param      rho_c  central density
* @param      n      polytropic index at the center of the star
* @param      pt     particle carrying the rest passive thermodynamic props
*/
std::pair<double, double>
lane_emden_RK4(const double m, const double s, const double th,
    const double dth, const double rho_c, const double n, body & pt) {
  double th12 = th + dth/2;
  double th2  = th + dth;
  auto [km_1, ks_1] = dms_dth(m,s,th,rho_c,n, pt);
  auto [km_2, ks_2] = dms_dth(m + dth/2*km_1,s + dth/2*ks_1,th12,rho_c,n,pt);
  auto [km_3, ks_3] = dms_dth(m + dth/2*km_2,s + dth/2*ks_2,th12,rho_c,n,pt);
  auto [km_4, ks_4] = dms_dth(m + dth*km_3,s + dth*ks_3,th2,rho_c,n,pt);

  double m_ret = m + dth/6*(km_1 + 2*km_2 + 2*km_3 + km_4);
  double s_ret = s + dth/6*(ks_1 + 2*ks_2 + 2*ks_3 + ks_4);

  return {m_ret, s_ret};
}

/**
* @brief      Lane-Emden solver
*
* @param      rho_c        central density
* @param      p_c          central pressure
* @param      Nr           radial resolution: number of points
* @param      rad_arr      radial grid output
* @param      rho_arr      density array
* @param      mass_arr     mass array
* @param      drhodr_arr   density derivative wrt r
*/
void
solve(const int Nr, std::vector<double> & rad_arr,
    std::vector<double> & rho_arr, std::vector<double> & mass_arr,
    std::vector<double> & drhodr_arr) {

  using namespace param;
  using namespace eos;
  const double rho_c = rho_initial;

  body pt0;
  pt0.setDensity(rho_initial);
  pt0.setPressure(pressure_initial);
  pt0.setAbar(initial_abar);
  pt0.setElectronfraction(initial_zbar/initial_abar);
  pt0.setTemperature(initial_temp);

  compute_internal_energy(pt0);
  compute_entropy(pt0);
  compute_pressure(pt0);
  compute_soundspeed(pt0);
  compute_internal_energy(pt0);
  const double p_c = pt0.getPressure();
  const double u_c = pt0.getInternalenergy();
  double cs = pt0.getSoundspeed();
  double CLIGHT2 = phys::clight * phys::clight;
  double dPdrho_c = cs*cs;
  if (lane_emden_isothermal) {
    if (get_dpdrho_at_temp == nullptr) {
      log_one(error) << "isothermal option not implemented for this EoS\n";
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
    dPdrho_c = get_dpdrho_at_temp(pt0);
  }

  // rho = rho_c * theta**n
  double gam = rho_c/p_c*dPdrho_c;
  double n = 1./(gam - 1.);

  // pseudo polytropic EOS for first step
  double K_c = p_c / pow(rho_c, gam);

  // useful constant for the first step
  double alpha = 4*M_PI*phys::GN / (K_c*(n+1)*pow(rho_c,(1.0/n)));

  // allocate arrays
  rad_arr.resize(Nr);
  rho_arr.resize(Nr);
  mass_arr.resize(Nr);
  drhodr_arr.resize(Nr);

  // start the solver
  // lane_emden_rho_atm the atmospheric pressure as the minimum pressure for integration
  // lane_emden_firststep can be any small number to prevent singularity (i.e. 1e-7)
  // double theta_min = lane_emden_firststep / (double)Nr;
  double theta_min = pow(lane_emden_rho_atm / rho_c, 1.0/n);
  double theta_step = - (1.0 - theta_min)/(double)(Nr - 1);
  std::vector<double> theta_arr(Nr);
  for(int i = 0; i < Nr; i++) {
      theta_arr[i] = 1.0 + i*theta_step;
  }

  // first step is approximated with polytropic EOS with const rho = rho_c, which gives
  // dm = 4*pi/3*rho_c*dr**3, ds = -6.0/(alpha*rho_c) * theta_step;
  double GR_cor_s_init = 1.0;
  double GR_cor_m_init = 1.0;
  if(tov_correction){
    GR_cor_s_init = 1.0/((1+(u_c*rho_c+p_c)/(CLIGHT2 * rho_c))*(1+(3*p_c)/(rho_c*CLIGHT2)));
    GR_cor_m_init = (1+u_c/CLIGHT2);
  }
  double s_init =-6./(alpha*rho_c) * theta_step * GR_cor_s_init;
  double m_init = 4.*M_PI/3*sqrt(CU(s_init)) * rho_c * GR_cor_m_init;
  double theta_cur = 1.;

  std::vector<double> s_arr(Nr);
  std::vector<double> m_arr(Nr);
  s_arr[0] = 0.;
  m_arr[0] = 0.;
  s_arr[1] = s_init;
  m_arr[1] = m_init;

  // RK4 for integration of the two ODEs
  for(int i = 1; i < Nr - 2; i++) {
    theta_cur = theta_arr[i];
    auto [first, second] = lane_emden_RK4(m_arr[i],s_arr[i],
        theta_arr[i], theta_step, rho_c, n, pt0);
    m_arr[i+1] = first;
    s_arr[i+1] = second;
  }

  // recursive refinement integration for the last step for singularity at theta = 0
  // the integration takes 20 iterations with the step size going half at each iteration
  double s_last = s_arr[Nr-2];
  double m_last = m_arr[Nr-2];
  double theta_last = theta_arr[Nr-2];
  double dtheta = theta_step;
  for(int i = 0; i < 20; i++){
    if(i < 19){
      dtheta *= 0.5;
    }
    auto [dm_dth_last, ds_dth_last] = dms_dth(m_last,s_last,theta_last,rho_c,n, pt0);
    s_last += dtheta*ds_dth_last;
    m_last += dtheta*dm_dth_last;
    theta_last += dtheta;
  }
  m_arr[Nr-1] = m_last;
  s_arr[Nr-1] = s_last;

  // Finally!
  double M_star = m_arr[Nr-1];
  double R_star = sqrt(s_arr[Nr-1]);

  // Output stellar parameters to log info
  log_one(info) << "\nLane-Emden solver:\n"
      << std::scientific << std::setprecision(12)
      << " - mass:    "<< M_star<< " [" << MASS_UNIT_STR << "]  = "
                       << M_star/phys::Msun << " [Msun]\n"
      << " - radius:  "<< R_star<< " [" << LENGTH_UNIT_STR << "] = "
                       << R_star/phys::Rsun << " [Rsun]\n"
      << " - central density:  " << rho_c << " [" << DENSITY_UNIT_STR << "]\n"
      << " - central pressure:  " << p_c << " [" << PRESSURE_UNIT_STR << "]\n"
      << std::endl;

  // RESETS param::sphere_radius to the value that has been found
  SET_PARAM(sphere_radius, R_star);

  // normalization constants
  const double
    rho_norm = M_star/CU(R_star),
    drhodr_norm = rho_norm/R_star;

  // normalize arrays to unit mass and unit radius
  for(int i = 0; i < Nr; i++){
    using namespace eos;
    double m = m_arr[i];
    double r = sqrt(s_arr[i]);
    double rho = rho_c * pow(theta_arr[i],n);
    pt0.setDensity(rho);
    compute_soundspeed(pt0);
    double cs = pt0.getSoundspeed();
    double dPdrho = lane_emden_isothermal
                  ? get_dpdrho_at_temp(pt0)
                  : cs*cs;
    double drhodr = -phys::GN*m*rho/(r*r * dPdrho);
    mass_arr[i] = m / M_star;
    rad_arr[i] = r / R_star;
    rho_arr[i] = rho / rho_norm;
    drhodr_arr[i] = (i ? drhodr/drhodr_norm : 0.);
  }

  rad_arr[0] = 0.;
  rho_arr[0] = rho_c / rho_norm;
  mass_arr[0] = 0.;
  drhodr_arr[0] = 0.;

  rad_arr[Nr - 1] = 1.;
  rho_arr[Nr - 1] = 0.;
  mass_arr[Nr - 1] = 1.;
  drhodr_arr[Nr - 1] = 0.;

  // Output the density profile with file name "lane_emden_output_profile" using cstdio
  // if string is empty (zero length), do not output profile
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if(rank == 0 and strlen(lane_emden_output_profile) > 0) {

    log_one(info) << "Generating output density profile in "
                  << lane_emden_output_profile << std::endl;

    // if the file already exists, issue a warning and overwrite it
    if(access(lane_emden_output_profile, F_OK ) != -1)
	  log_one(warn) << "File exists: overwriting "
                    << lane_emden_output_profile << std::endl;

    // create header
    std::ostringstream oss_header;
    oss_header << "# Stellar parameters:\n" << std::setprecision(12)
      << "#  - mass:    " << M_star << " [" << MASS_UNIT_STR <<"] = "
                          << (M_star/phys::Msun) << " [Msun]\n"
      << "#  - radius:  " << R_star << " [" << LENGTH_UNIT_STR << "] = "
                          << (R_star/phys::Rsun) << " [Rsun]\n"
      << "#  - central density:   " << rho_c << " [" << DENSITY_UNIT_STR << "]\n"
      << "#  - central pressure:  " << p_c << " [" << PRESSURE_UNIT_STR << "]\n"
      << "#  - electron fraction: " << (initial_zbar/initial_abar) << "\n"
      << "#\n"
      << "# Equation of state: " << eos_type_decode[(int)eos_type]
      << std::endl;

    switch (eos_type) {
      case param::eos_ideal:
      case param::eos_polytropic:
        oss_header << "#  - poly_gamma = " << poly_gamma << std::endl;
        break;

      case param::eos_ppt:
        oss_header << "#  - poly_gamma  = " << poly_gamma << "\n"
                   << "#  - poly_gamma2 = " << poly_gamma2 << "\n"
                   << "#  - ppt_density_thr = " << ppt_density_thr
                   << std::endl;
        break;

      case param::eos_helmholtz:
      case param::eos_stellar_collapse:
        oss_header << "#  - eos_tab_file_path: "
                   << "\"" << eos_tab_file_path << "\""
                   << std::endl;
    }

    std::ofstream out(lane_emden_output_profile);
    out << oss_header.str();

    // output the profile data using format below
    for(int i = 0; i < Nr; i++) {
      out << std::scientific << std::setprecision(12)
          << rad_arr[i]  << " " << rho_arr[i] << " "
          << mass_arr[i] << " " << drhodr_arr[i]
          << "\n";
    }
    out << std::flush;
    out.close();

    // check that the file has been written; if not: complain and exit
    if(access(lane_emden_output_profile, F_OK ) == -1) {
      log_one(error) << "\n Density profile cannot be created" << std::endl;
      MPI_Abort (MPI_COMM_WORLD, -1);
    }
  }
} // solve(..)

} // namespace lane_emden


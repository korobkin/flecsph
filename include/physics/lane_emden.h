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

#include "eos.h"
#include "body.h"

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
  double rho = rho_c * pow(th, n);
  pt.setDensity(rho);
  eos::compute_soundspeed(pt);
  double cs = pt.getSoundspeed();
  double dPdrho_S = cs*cs;
  double dsdth = -2*n*sqrt(s*s*s)/(GNEWT * m * th) * dPdrho_S;
  double dmdth = dsdth * 2*M_PI*sqrt(s)*rho;
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
  const double rho_c = rho_initial;
  const double p_c = pressure_initial;

  body pt0;
  pt0.setDensity(rho_initial);
  pt0.setPressure(pressure_initial);
  pt0.setAbar(initial_abar);
  pt0.setElectronfraction(initial_zbar/initial_abar);

  eos::compute_entropy(pt0);
  eos::compute_soundspeed(pt0);
  double cs = pt0.getSoundspeed();
  double dPdrho_c = cs*cs;

  // rho = rho_c * theta**n
  double gam = rho_c/p_c*dPdrho_c;
  double n = 1./(gam - 1.);

  // pseudo polytropic EOS for first step
  double K_c = p_c / pow(rho_c, gam);

  // useful constant for the first step
  double alpha = 4*M_PI*GNEWT / (K_c*(n+1)*pow(rho_c,(1.0/n)));

  // allocate arrays
  rad_arr.resize(Nr);
  rho_arr.resize(Nr);
  mass_arr.resize(Nr);
  drhodr_arr.resize(Nr);

  //start the solver
  double theta_min = 1e-7 / (double)Nr;  // TODO: make 1e-7 into a parameter
  double theta_step = - (1.0 - theta_min)/(double)(Nr - 1);
  std::vector<double> theta_arr(Nr);
  for(int i = 0; i < Nr; i++) {
      theta_arr[i] = 1.0 + i*theta_step;
  }

  // first step is approximated with polytropic EOS with const rho = rho_c, which gives
  // dm = 4*pi/3*rho_c*dr**3, ds = -6.0/(alpha*rho_c) * theta_step;
  double s_init =-6./(alpha*rho_c) * theta_step;
  double m_init = 4.*M_PI/3*sqrt(CU(s_init)) * rho_c;
  double theta_cur = 1.;

  std::vector<double> s_arr(Nr);
  std::vector<double> m_arr(Nr);
  s_arr[0] = 0.;
  m_arr[0] = 0.;
  s_arr[1] = s_init;
  m_arr[1] = m_init;

  // RK4 for integration of the two ODEs
  for(int i = 1; i < Nr - 1; i++) {
    theta_cur = theta_arr[i];
    auto [first, second] = lane_emden_RK4(m_arr[i],s_arr[i],
        theta_arr[i], theta_step, rho_c, n, pt0);
    m_arr[i+1] = first;
    s_arr[i+1] = second;
  }

  // Finally!
  double M_star = m_arr[Nr-1];
  double R_star = sqrt(s_arr[Nr-1]);

  // Output stellar parameters to log info
  log_one(info) << "\nLane-Emden solver:\n"
      << std::scientific << std::setprecision(12)
      << " - mass:    "<< M_star<< " [g]  = "<<M_star/M_SUN_CGS<< " [Msun]\n"
      << " - radius:  "<< R_star<< " [cm] = "<<R_star/R_SUN_CGS<< " [Rsun]\n"
      << " - central density:  " << rho_c << " [g/cm^3]\n"
      << " - central pressure:  " << p_c << " [dynes/cm^2]\n"
      << std::endl;

  // RESETS param::sphere_radius to the value that has been found
  SET_PARAM(sphere_radius, R_star);

  // normalization constants
  const double
    rho_norm = M_star/CU(R_star),
    drhodr_norm = rho_norm/R_star;

  // normalize arrays to unit mass and unit radius
  for(int i = 0; i < Nr; i++){
    double m = m_arr[i];
    double r = sqrt(s_arr[i]);
    double rho = rho_c * pow(theta_arr[i],n);
    pt0.setDensity(rho);
    eos::compute_soundspeed(pt0);
    double cs = pt0.getSoundspeed();
    double dPdrho_S = cs*cs;
    double drhodr = -GNEWT*m*rho/(r*r * dPdrho_S);
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

  // // UNCOMMENT for quick-and-dirty profile output to stdout
  // // TODO: 1. add a parameter: string 'lane_emden_output_profile'
  // //       2. only output from MPI rank 0
  // //       2. if string is empty (zero length), do not output profile;
  // //       3. if string is non-empty, assume it contains profile file name;
  // //       4. attempt to create file with that name;
  // //       5. if the file already exists, issue a warning and overwrite it;
  // //       6. check that the file has been successfully created;
  // //       7. output the header (make sure to correctly specify EOS)
  // 
  // printf ("# Stellar parameters:\n");
  // printf ("#  - mass:    %12.12e [g]\n", M_star);
  // printf ("#  - radius:  %12.12e [cm]\n", R_star);
  // printf ("#  - central density:  %12.12e [g/cm^3]\n", rho_c);
  // printf ("#  - central pressure:  %12.12e [dynes/cm^2]\n", p_c);
  // printf ("#\n");
  // printf ("# Equation of state: zero-temperature WD\n");
  // 
  // //       8. output the profile data using format below
  // for(int i = 0; i < Nr; i++){
  //   printf("%19.12e %19.12e %19.12e %19.12e\n",
  //       rad_arr[i], rho_arr[i], mass_arr[i], drhodr_arr[i]);
  // }

} // solve(..)

} // namespace lane_emden


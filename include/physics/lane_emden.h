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
#include "eos.h"
#include "body.h"
#include "params.h"
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
  eos::compute_pressure(pt);
  //eos::compute_soundspeed(pt);
  eos::compute_internal_energy(pt);
  double p = pt.getPressure();
  double u = pt.getInternalenergy();
  //double cs = pt.getSoundspeed();
  const double CLIGHT2 = C_LIGHT_CGS * C_LIGHT_CGS;
  double dPdrho_S = eos::get_dPdrho(pt);
  // tov correction terms
  double GR_cor_ds = 1.0;
  double GR_cor_dm = 1.0;
  if(param::tov_correction){
    double GR_cor_ds1 = (1 + (rho*u + p) / (CLIGHT2*rho));
    double GR_cor_ds2 = (1 + (4*M_PI*sqrt(s*s*s)*p)/(m*CLIGHT2));
    double GR_cor_ds3 = (1 - (2*GNEWT*m)/(sqrt(s)*CLIGHT2));
    
    GR_cor_ds = GR_cor_ds3 / (GR_cor_ds1 * GR_cor_ds2);
    GR_cor_dm = (1 + u/CLIGHT2);
  }

  double dsdth = -2*n*sqrt(s*s*s)/(GNEWT * m * th) * dPdrho_S * GR_cor_ds;
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
solve(const int Nr, std::vector<double> & rad_arr, std::vector<double> & rho_arr,
                    std::vector<double> & mass_arr, std::vector<double> & drhodr_arr,
                    std::vector<double> & alpha2_arr, std::vector<double> & dalpha2dr_arr,
                    std::vector<double> & beta2_arr, std::vector<double> & dbeta2dr_arr) {

  using namespace param;
  const double rho_c = rho_initial;

  body pt0;
  pt0.setDensity(rho_initial);
  pt0.setPressure(pressure_initial);
  pt0.setAbar(initial_abar);
  pt0.setElectronfraction(initial_zbar/initial_abar);
  pt0.setTemperature(initial_temp);

  eos::compute_internal_energy(pt0);
  eos::compute_entropy(pt0);
  eos::compute_pressure(pt0);
  eos::compute_soundspeed(pt0);
  eos::compute_internal_energy(pt0);
  const double p_c = pt0.getPressure();
  const double cs_c = pt0.getSoundspeed();
  const double u_c = pt0.getInternalenergy();
  double CLIGHT2 = C_LIGHT_CGS * C_LIGHT_CGS;
  double dPdrho_c = eos::get_dPdrho(pt0);

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
  alpha2_arr.resize(Nr);
  dalpha2dr_arr.resize(Nr);
  beta2_arr.resize(Nr);
  dbeta2dr_arr.resize(Nr);
  // start the solver
  // lane_emden_rho_atm the atmospheric density as the minimum density for integration (at theta = theta_min)
  double theta_min = pow(lane_emden_rho_atm / rho_c, 1.0/n);
  
  // theta_cutoff roughly seperates the integration to the crust and core
  double theta_cutoff = 0.001;
  
  // the core has 0.9*Nr regularly-spaced grid, the crust has 0.1*Nr logarithmicly-spaced grids
  int Nr_core = (int)(0.9*Nr);
  int Nr_crust = Nr - Nr_core;
  double theta_core_step = - (1.0 - theta_cutoff)/(double)(Nr_core - 1);
  double theta_crust_step_log = pow((theta_min / theta_cutoff),1.0/(Nr_crust-1));
  std::vector<double> theta_arr(Nr);
  for(int i = 0; i < Nr; i++) {
      if(i < Nr_core){
        theta_arr[i] = 1.0 + i*theta_core_step;
      } else {
        theta_arr[i] = theta_arr[i-1] * theta_crust_step_log;
      }
  }
  double theta_step = theta_arr[1] - theta_arr[0];

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
  for(int i = 1; i < Nr - 1; i++) {
    theta_cur = theta_arr[i];
    theta_step = theta_arr[i+1] - theta_arr[i];
    auto [first, second] = lane_emden_RK4(m_arr[i],s_arr[i],
        theta_arr[i], theta_step, rho_c, n, pt0);
    m_arr[i+1] = first;
    s_arr[i+1] = second;
  }

  // recursive refinement integration for the last step for singularity at theta = 0
  // the integration takes 20 iterations with the step size going half at each iteration
//  double s_last = s_arr[Nr-2];
//  double m_last = m_arr[Nr-2];
//  double theta_last = theta_arr[Nr-2];
//  double dtheta = theta_step;
//  for(int i = 0; i < 20; i++){
//    if(i < 19){
//      dtheta *= 0.5;
//    }
//    auto [dm_dth_last, ds_dth_last] = dms_dth(m_last,s_last,theta_last,rho_c,n, pt0);
//    s_last += dtheta*ds_dth_last;
//    m_last += dtheta*dm_dth_last;
//    theta_last += dtheta;
//  }
//  m_arr[Nr-1] = m_last;
//  s_arr[Nr-1] = s_last;

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
      << " - central soundspeed:  " << cs_c << " [cm/s]\n"
      << std::endl;

  // RESETS param::sphere_radius to the value that has been found
  SET_PARAM(sphere_radius, R_star);
  SET_PARAM(sphere_mass, M_star);
  // fixedGR background metric required alpha^2, beta^2, d(alpha^2)/dr, d(beta^2)/dr
  // -gtt = alpha^2 = exp^(2Phi), grr = beta^2 = exp^(2Lambda)
  // beta^2 =  1.0/(1-2Gm/rc^2), dPhi/dr = G/c^2((m+4pi*r^3p/c^2)/(r*(r-2Gm/c^2))
  
  //std::vector<double> beta2(Nr);
  //std::vector<double> dbeta2dr(Nr);
  std::vector<double> Lambda(Nr);
  //std::vector<double> alpha2(Nr);
  //std::vector<double> dalpha2dr(Nr);
  std::vector<double> Phi(Nr);
  std::vector<double> dPhidr(Nr);
  std::vector<double> r_arr(Nr);
  for(int i = 0; i < Nr; i++){
    r_arr[i] = sqrt(s_arr[i]);
  }
  
  Phi[Nr-1] = log(1-2*GNEWT*M_star/(R_star*CLIGHT2))/2.0;
  for(int i = 0; i < Nr; i++){
    pt0.setDensity(rho_arr[i]);
    eos::compute_pressure(pt0);
    //eos::compute_soundspeed(pt0);
    eos::compute_internal_energy(pt0);
    double p_cur = pt0.getPressure();
    double u_cur = pt0.getInternalenergy();
    // double cs_cur = pt0.getSoundspeed();
    double rsh = 2*GNEWT*m_arr[i]/CLIGHT2;

    beta2_arr[i] = 1.0/(1 - rsh/r_arr[i]);
    double dmdr = 4*M_PI*pow(r_arr[i],2)*rho_arr[i]*(1+u_cur/CLIGHT2);
    dbeta2dr_arr[i] = rsh/r_arr[i] * pow(1-rsh/r_arr[i],-2.0)*(dmdr/m_arr[i] - 1.0/r_arr[i]);
    Lambda[i] = log(beta2_arr[i])/2;

    dPhidr[i] = rsh * (1 + 4*M_PI*pow(r_arr[i],3)*p_cur/(m_arr[i] * CLIGHT2)) / (r_arr[i]*(r_arr[i] - rsh));
  }
  // integrate Phi using traoezoidal rule from surface back to the core
  for(int i = Nr-2; i >= 0; i--){
    Phi[i] = Phi[i+1] - (r_arr[i+1] - r_arr[i]) * (dPhidr[i] + dPhidr[i+1])/2;
    alpha2_arr[i] = exp(2*Phi[i]);
    dalpha2dr_arr[i] = 2*alpha2_arr[i]*dPhidr[i];
  }

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
    //eos::compute_soundspeed(pt0);
    eos::compute_pressure(pt0);
    eos::compute_internal_energy(pt0);
    // double cs = pt0.getSoundspeed();
    double p = pt0.getPressure();
    double u = pt0.getInternalenergy();
    double dPdrho_S = eos::get_dPdrho(pt0);
    double drhodr = -GNEWT*m*rho/(r*r * dPdrho_S);
    mass_arr[i] = m / M_star;
    rad_arr[i] = r / R_star;
    rho_arr[i] = rho / rho_norm;
    //p_arr[i] = p;
    //u_arr[i] = u;
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
      << "#  - mass:    " << M_star << " [g] = "
                          << (M_star/M_SUN_CGS) << " [Msun]\n"
      << "#  - radius:  " << R_star << " [cm] = "
                          << (R_star/R_SUN_CGS) << " [Rsun]\n"
      << "#  - central density:   " << rho_c << " [g/cm^3]\n"
      << "#  - central pressure:  " << p_c << " [dynes/cm^2]\n"
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


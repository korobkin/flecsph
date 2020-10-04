/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2018 Triad National Security, LLC
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
 * @file eos.h
 * @brief Namespace for both analytic and tabulated EOS
 *        for pressure computation and sound speed
 */

#pragma once

#include <vector>

#include "params.h"
#include "tree.h"
#include "utils.h"
#include <boost/algorithm/string.hpp>

#ifdef ENABLE_DEBUG_EOS
#  define _DEBUG_EOS_
#  warning "Debug mode for equations of state"
#endif

#define SQ(x) ((x) * (x))
#define CU(x) ((x) * (x) * (x))
#define QU(x) ((x) * (x) * (x) * (x))
namespace eos {
using namespace param;

template<>
class eos_t<param::eos_ppt>{
  static double rho_thr;    // density threshold1
  static double rho_thr2;    // density threshold2
public:
  /**
  * @brief      Compute adiabatic invariant (a function of entropy)
  *             from density and pressure.
  *             In the piecewise-polytropic EOS, we pick adiabatic
  *             invariant to be the constant over the first polytropic
  *             segment K1:
  *
  *              P(rho) = K1\rho^\Gamma1 + K2\rho^\Gamma2
  *
  * @param      particle
  */
  static void
  compute_entropy(body & particle){
    const double rho = particle.getDensity(),
                 P   = particle.getPressure();
    double K1 = 0.0;
    if (ppt_pressure_thr > 0) {
      K1 = ppt_pressure_thr/pow(rho_thr, poly_gamma);
    }
    else if (rho < rho_thr) {
      K1 = P/pow(rho, poly_gamma);
    }
    else if (rho < rho_thr2) {
      double K2 = P/pow(rho, poly_gamma2);
      K1 = K2*pow(rho_thr , poly_gamma2 - poly_gamma);
    } else {
      double K3 = P/pow(rho, poly_gamma3);
      double K2 = K3*pow(rho_thr2, poly_gamma3 - poly_gamma2);
      K1 = K2*pow(rho_thr , poly_gamma2 - poly_gamma );
    }
    particle.setEntropy(K1);
  }

  /**
  * @brief      Initialize equation of state:
  *             - it ppt_pressure_thr is specified, reset pressure_initial
  */
  static void init() {
    eos_t<param::eos_ppt>::rho_thr  = param::ppt_density_thr;
    eos_t<param::eos_ppt>::rho_thr2 = param::ppt_density_thr2;
    if (ppt_pressure_thr > 0) {
      double K1 = ppt_pressure_thr/pow(ppt_density_thr, poly_gamma);
      body pt;
      pt.setDensity(rho_initial);
      pt.setEntropy(K1);
      compute_pressure(pt);
      SET_PARAM(pressure_initial, pt.getPressure());
    }
  }

  /**
  * @brief      Compute the pressure for piecewise-polytrope EOS
  * @param      particle
  */
  static void
  compute_pressure(body & particle) {
    const double rho = particle.getDensity(),
                 K1  = particle.getEntropy();
    double P = 0.0;
    if (rho < rho_thr) {
      P = K1*pow(rho, poly_gamma);
    }
    else if (rho < rho_thr2) {
      double K2 = K1*pow(rho_thr , poly_gamma  - poly_gamma2);
      P = K2*pow(rho, poly_gamma2);
    } 
    else {
      double K2 = K1*pow(rho_thr , poly_gamma  - poly_gamma2);
      double K3 = K2*pow(rho_thr2, poly_gamma2 - poly_gamma3);
      P = K3*pow(rho, poly_gamma3);
    }
    particle.setPressure(P);
  }

  /**
  * @brief      Compute sound speed for piecewise polytropic eos
  *
  * @param      particle
  */
  static void
  compute_soundspeed(body & particle) {
    const double rho = particle.getDensity(),
                 K1  = particle.getEntropy(),
                 gam = (rho < rho_thr ? poly_gamma : (rho < rho_thr2 ? poly_gamma2 : poly_gamma3));
    double soundspeed = 0.;
    if (rho < rho_thr) {
      soundspeed = sqrt(K1*poly_gamma*pow(rho,poly_gamma - 1.));
    }
    else if (rho < rho_thr2) {
      double K2 = K1*pow(rho_thr, poly_gamma - poly_gamma2);
      soundspeed = sqrt(K2*poly_gamma2*pow(rho,poly_gamma2 - 1.));
    } else {
      double K2 = K1*pow(rho_thr , poly_gamma  - poly_gamma2);
      double K3 = K2*pow(rho_thr2, poly_gamma2 - poly_gamma3);
      soundspeed = sqrt(K3*poly_gamma3*pow(rho,poly_gamma3 - 1.));
    }
    particle.setSoundspeed(soundspeed);
  }

  /**
  * @brief      Empty function because EOS is temperature-agnostic
  *             Can be tied to internal energy via <A> and IG equation
  *
  * @param      particle
  */
  static void
  compute_temperature(body&){}

  /**
  * @brief      Compute specific internal energy
  *             Uses adiabatic invariant and density
  *
  * @param      particle
  */
  static void
  compute_internal_energy(body & particle) {
    const double rho = particle.getDensity(),
                 K1  = particle.getEntropy();
    double eps = 0.;
    if (rho < rho_thr) {
      eps = K1*pow(rho, poly_gamma - 1.)/(poly_gamma - 1.);
    }
    else if (rho < rho_thr2) {
      double K2 = K1*pow(rho_thr, poly_gamma - poly_gamma2);
      eps = K2*pow(rho,     poly_gamma2 - 1.)/(poly_gamma2 - 1.)
          - K2*pow(rho_thr, poly_gamma2 - 1.)/(poly_gamma2 - 1.)
          + K1*pow(rho_thr, poly_gamma  - 1.)/(poly_gamma  - 1.);
    } else {
      double K2 = K1*pow(rho_thr , poly_gamma  - poly_gamma2);
      double K3 = K2*pow(rho_thr2, poly_gamma2 - poly_gamma3);
      eps = K3*pow(rho,      poly_gamma3 - 1.)/(poly_gamma3 - 1.)
          - K3*pow(rho_thr2, poly_gamma3 - 1.)/(poly_gamma3 - 1.)
          + K2*pow(rho_thr2, poly_gamma2 - 1.)/(poly_gamma2 - 1.)
          - K2*pow(rho_thr , poly_gamma2 - 1.)/(poly_gamma2 - 1.)
          + K1*pow(rho_thr , poly_gamma  - 1.)/(poly_gamma  - 1.);
    }
    particle.setInternalenergy(eps);
  }

  compute_quantity_t compute_spct_given_rho_u = nullptr;

}; // class eos_t<param::eos_ppt>

// declare static member of a templated class
template<>
double eos_t<param::eos_ppt>::rho_thr;

template<>
double eos_t<param::eos_ppt>::rho_thr2;

} // namespace eos

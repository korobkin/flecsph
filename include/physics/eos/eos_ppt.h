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
  static constexpr int max_num_segments = 12;
  static int num_segments;
  static double gammas[max_num_segments];
  static double rho_thr[max_num_segments];  // density thresholds

public:
  /**
  * @brief      Initialize equation of state:
  *             - it ppt_pressure_thr is specified, reset pressure_initial
  */
  static void init() {
    // by default, have three segments
    num_segments = 3;
    rho_thr[0] = ppt_density_thr;
    rho_thr[1] = ppt_density_thr2;
    gammas[0] = poly_gamma;
    gammas[1] = poly_gamma2;
    gammas[2] = poly_gamma3;

    // Parameter ppt_pressure_thr corresponds to the pressure at first
    // density threshold; if it is specified, then all constants K1, K2, ..
    // are fixed. 
    // The followin segment recomputes pressure_initial at rho_initial:
    // pressure_initial = P(rho_initial)
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
  * @brief      Compute adiabatic invariant (a function of entropy)
  *             from density and pressure.
  *             In the piecewise-polytropic EOS, we pick adiabatic
  *             invariant to be the constant over the first polytropic
  *             segment K1:
  *
  *              P(rho) = K1\rho^\Gamma1
  *
  *             If ppt_pressure_thr is set, then K1 is computed using
  *             ppt_pressure_thr, ppt_density_thr and poly_gamma
  * @param      particle
  */
  static void
  compute_entropy(body & particle){
    const double rho = particle.getDensity(),
                 P   = particle.getPressure();
    double K1 = 0.0;
    if (ppt_pressure_thr > 0) {
      K1 = ppt_pressure_thr/pow(ppt_density_thr, poly_gamma);
    }
    else {
      for (int i = 0; i < num_segments; ++i) {
        if (rho < rho_thr[i] or i == num_segments - 1) {
          double K1 = P/pow(rho, gammas[i]);
          for (int j = i; j > 0; --j) 
            K1 *= pow(rho_thr[j - 1] , gammas[j] - gammas[j - 1]);
          break;
        }
      }
    }
    particle.setEntropy(K1);
  }

  /**
  * @brief      Compute the pressure for piecewise-polytrope EOS
  *             Uses density rho and entropy function K1
  * @param      particle
  */
  static void
  compute_pressure(body & particle) {
    double rho = particle.getDensity(),
           Kn  = particle.getEntropy();
    int i = 0;
    for (; i < num_segments - 1; ++i) {
      if (rho < rho_thr[i])
        break;
      else 
        Kn *= pow(rho_thr[i], gammas[i] - gammas[i + 1]);
    }
    particle.setPressure(Kn*pow(rho, gammas[i]));
  }

  /**
  * @brief      Compute sound speed for piecewise polytropic eos
  *             Uses density rho and entropy function K1
  *
  * @param      particle
  */
  static void
  compute_soundspeed(body & particle) {
    double rho = particle.getDensity(),
           Kn  = particle.getEntropy();
    int i = 0;
    for (; i < num_segments - 1; ++i) {
      if (rho < rho_thr[i])
        break;
      else 
        Kn *= pow(rho_thr[i], gammas[i] - gammas[i + 1]);
    }
    double cs = sqrt(Kn*gammas[i]*pow(rho, gammas[i] - 1.));
    particle.setSoundspeed(cs);
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
    double rho = particle.getDensity(),
           Kn  = particle.getEntropy(),
           delta_eps = 0.;
    int i = 0;
    for (; i < num_segments - 1; ++i) {
      if (rho < rho_thr[i])
        break;
      else {
        delta_eps += Kn*pow(rho_thr[i], gammas[i]  - 1.)/(gammas[i]  - 1.);
        Kn *= pow(rho_thr[i], gammas[i] - gammas[i + 1]);
        delta_eps -= Kn*pow(rho_thr[i], gammas[i+1]- 1.)/(gammas[i+1]- 1.);
      }
    }

    double eps = Kn*pow(rho, gammas[i] - 1.)/(gammas[i] - 1.) + delta_eps;
    particle.setInternalenergy(eps);
  }

  compute_quantity_t compute_spct_given_rho_u = nullptr;

}; // class eos_t<param::eos_ppt>

// declare static members of a templated class
template<>
int eos_t<eos_ppt>::num_segments;

template<>
double eos_t<eos_ppt>::gammas[eos_t<eos_ppt>::max_num_segments];

template<>
double eos_t<eos_ppt>::rho_thr[eos_t<eos_ppt>::max_num_segments];

} // namespace eos

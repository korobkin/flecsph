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

    // ppt fits to various cold NS EoS from Read'09 (arXiv:0812.2163)
    if (boost::equals(ppt_eos_fit, "crust")) {
      // crust fit from Table II in the Appendix C
      num_segments = 4;
      rho_thr[0] = 2.44034e+07;
      rho_thr[1] = 3.78358e+11;
      rho_thr[2] = 2.62780e+12;

      gammas[0] = 1.58425;
      gammas[1] = 1.28733;
      gammas[2] = 0.62223;
      gammas[3] = 1.35692;

      SET_PARAM(ppt_pressure_thr, 3436.745395);
      SET_PARAM(poly_gamma, gammas[0]);

    }
    else if (not boost::equals(ppt_eos_fit, "none")) {
      struct ppt_fit {
        char name[8];
        double lg_p1;
        double gamma1, gamma2, gamma3;
      };

      // fits from Table III in the Appendix C
      ppt_fit fits[] = {
        {  "PAL6",  34.380,  2.227,  2.189,  2.159},
        {   "SLy",  34.384,  3.005,  2.988,  2.851},
        {   "AP1",  33.943,  2.442,  3.256,  2.908},
        {   "AP2",  34.126,  2.643,  3.014,  2.945},
        {   "AP3",  34.392,  3.166,  3.573,  3.281},
        {   "AP4",  34.269,  2.830,  3.445,  3.348},
        {   "FPS",  34.283,  2.985,  2.863,  2.600},
        {  "WFF1",  34.031,  2.519,  3.791,  3.660},
        {  "WFF2",  34.233,  2.888,  3.475,  3.517},
        {  "WFF3",  34.283,  3.329,  2.952,  2.589},
        {  "BBB2",  34.331,  3.418,  2.835,  2.832},
        {"BPAL12",  34.358,  2.209,  2.201,  2.176},
        {   "ENG",  34.437,  3.514,  3.130,  3.168},
        {  "MPA1",  34.495,  3.446,  3.572,  2.887},
        {   "MS1",  34.858,  3.224,  3.033,  1.325},
        {   "MS2",  34.605,  2.447,  2.184,  1.855},
        {  "MS1b",  34.855,  3.456,  3.011,  1.425},
        {    "PS",  34.671,  2.216,  1.640,  2.365},
        {   "GS1",  34.504,  2.350,  1.267,  2.421},
        {   "GS2",  34.642,  2.519,  1.571,  2.314},
        {"BGN1H1",  34.623,  3.258,  1.472,  2.464},
        {  "GNH3",  34.648,  2.664,  2.194,  2.304},
        {    "H1",  34.564,  2.595,  1.845,  1.897},
        {    "H2",  34.617,  2.775,  1.855,  1.858},
        {    "H3",  34.646,  2.787,  1.951,  1.901},
        {    "H4",  34.669,  2.909,  2.246,  2.144},
        {    "H5",  34.609,  2.793,  1.974,  1.915},
        {    "H6",  34.593,  2.637,  2.121,  2.064},
        {    "H7",  34.559,  2.621,  2.048,  2.006},
        {  "PCL2",  34.507,  2.554,  1.880,  1.977},
        {  "ALF1",  34.055,  2.013,  3.389,  2.033},
        {  "ALF2",  34.616,  4.070,  2.411,  1.890},
        {  "ALF3",  34.283,  2.883,  2.653,  1.952},
        {  "ALF4",  34.314,  3.009,  3.438,  1.803}
      };

      num_segments = 3;
      rho_thr[0] = exp10(14.7);
      rho_thr[1] = 1e+15;
      const size_t num_fits = sizeof(fits)/sizeof(ppt_fit);
      int i = 0;
      for (; i < num_fits; ++i) {
        if (boost::equals(ppt_eos_fit, fits[i].name)) {
          log_one(info) << "setting piecewise-polytropic fit values for EoS \""
                        << fits[i].name << "\"" << std::endl;
          gammas[0] = fits[i].gamma1;
          gammas[1] = fits[i].gamma2;
          gammas[2] = fits[i].gamma3;
          SET_PARAM(poly_gamma, gammas[0]);
          SET_PARAM(ppt_pressure_thr, exp10(fits[i].lg_p1));
          break;
        }
      }
      if (i == num_fits) { // not found
        log_one(error) << "ppt_eos_fit: unknown value \""
                       <<  ppt_eos_fit << "\"" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
      }
    }

    // Parameter ppt_pressure_thr corresponds to the pressure at first
    // density threshold; if it is specified, then all constants K1, K2, ..
    // are fixed.
    // The followin segment recomputes pressure_initial at rho_initial:
    // pressure_initial = P(rho_initial)
    if (ppt_pressure_thr > 0) {
      double K1 = ppt_pressure_thr/pow(rho_thr[0], gammas[0]);
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

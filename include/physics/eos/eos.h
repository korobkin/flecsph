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

#include "eos_utils.h"
#include "eos_consts.h"
#include "eos_stellar_collapse.h"
#include "eos_helm.h"

namespace eos {
using namespace param;

constexpr double square(const double& x){
  return ((x) * (x));
}
constexpr double cube(const double& x){
  return ((x) * (x) * (x));
}
constexpr double quartic(const double& x){
  return ((x) * (x) * (x) * (x));
}


template<>
class eos_t<param::eos_ideal>{

public:
  static void init(body & particle){}

  static void read_data(){}

  /**
  * @brief      Computes pressure using the density and internal energy
  *
  * @param      particle
  */
  static void compute_pressure(body& particle){
    double pressure =
      (poly_gamma - 1.0) * particle.getDensity() * particle.getInternalenergy();
    particle.setPressure(pressure);
  }

  /**
  * @brief      Sound speed from specific internal energy
  *
  * @param      particle
  */
  static void compute_soundspeed(body & particle) {
    const double eps = particle.getInternalenergy();
    double soundspeed = sqrt(poly_gamma*(poly_gamma - 1)*eps);
    particle.setSoundspeed(soundspeed);
  }

  /**
  * @brief      Compute temperature via ideal gas formula
  *             Uses abar and specific internal energy
  *
  * @param      particle
  */
  static void compute_temperature(body & particle) {
    const double abar = particle.getAbar(),
                 eps  = particle.getInternalenergy();
    double T = AMU/KBOL*abar*(poly_gamma - 1)*eps;
    particle.setTemperature(T);
  }

  /**
  * @brief      Compute specific internal energy
  *             Uses adiabatic invariant and density
  *
  * @param      particle
  */
  static void
  compute_internal_energy(body & particle) {
    const double rho = particle.getDensity(),
                 K   = particle.getAdiabatic();
    double eps = K*pow(rho, poly_gamma - 1)/(poly_gamma - 1);
    particle.setInternalenergy(eps);
  }

  /**
  * @brief      Compute specific internal energy
  *             Uses abar and temperature
  *
  * @param      particle
  */
  static void
  compute_internal_energy_from_temperature(body & particle) {
    const double abar = particle.getAbar(),
                 T    = particle.getTemperature();
    double eps = T*KBOL/(AMU*abar*(poly_gamma - 1));
    particle.setInternalenergy(eps);
  }

}; // ...<eos_ideal>

template<>
class eos_t<param::eos_polytropic>{

public:
  /**
  * @brief      Initialized adiabatic invariant from initial conditions
  *
  * @param      particle
  */
  static void init(body & particle){
    const double rho = particle.getDensity(),
                 P = particle.getPressure();
    double K = P/pow(rho,poly_gamma);
    particle.setAdiabatic(K);
    return;
  }

  static void read_data(){}

  /**
  * @brief      Compute pressure from density using polytrope
  *             P(\rho) = A*\rho^\Gamma
  *
  * @param      particle
  */
  static void compute_pressure(body & particle) {
    const double rho = particle.getDensity(),
                 K = particle.getAdiabatic();
    particle.setPressure(K*pow(rho, poly_gamma));
  }

  /**
  * @brief      Compute sound speed for ideal fluid or polytropic eos
  *             cs = sqrt{ \Gamma\rho^(\Gamma-2) }
  *
  * @param      particle
  */
  static void compute_soundspeed(body & particle) {
    const double rho = particle.getDensity(),
                 K = particle.getAdiabatic();
    double soundspeed = sqrt(K*poly_gamma*pow(rho, poly_gamma - 2));
    particle.setSoundspeed(soundspeed);
  }

  /**
  * @brief      For polytropic equation of state, the temperature is
  *             decoupled from density or pressure, so this function does
  *             nothing
  *
  * @param      particle
  */
  static void
  compute_temperature(body& particle){}

  /**
  * @brief      Compute specific internal energy
  *             Uses adiabatic invariant and density
  *
  * @param      particle
  */
  static void
  compute_internal_energy(body & particle) {
    const double rho = particle.getDensity(),
                 K   = particle.getAdiabatic();
    double eps = K*pow(rho, poly_gamma - 1)/(poly_gamma - 1);
    particle.setInternalenergy(eps);
  }

}; // ...<eos_polytropic>


/**
* @brief      Equation of state for a cold white dwarf.
*             The pressure function psi(x)
*
*               psi(x) = (x*(2*x^2 - 3) * sqrt(1 + x^2) + 3*asinh(x))
*
*             can be fit reasonably well with a piecewise polytrope:
*                         | A1 x^5, if x < x0
*               psi(x) = <
*                         | A2 x^4, if x > x0
*             where x0 = 1.25, A1 = 1.6 and A2 = 2.0.
*
*/
template<>
class eos_t<param::eos_wd>{

  // pressure function constants
  static constexpr double A_wd = 6.00288e22;
  static constexpr double B_wd_nm = 9.81011e5;

  // constants of the piecewoise-polytrope fit to the pressure function
  static constexpr double ppt_x0 = 1.25;
  static constexpr double ppt_A1 = 1.6;
  static constexpr double ppt_A2 = 2.0;

public:
  static void init(body & particle){
    compute_internal_energy(particle);
  }

  static void read_data(){}

  static inline double
  pressure_from_rhoYe(double rho, double Ye) {
    double x = cbrt(rho*Ye/B_wd_nm);
    double x2 = x*x;
    return A_wd*(x*(2*x2 - 3) * sqrt(x2 + 1) + 3*asinh(x));
  }

  static inline double
  soundspeed_from_rhoYe(double rho, double Ye) {
    double x = cbrt(rho*Ye/B_wd_nm);
    double x2 = x*x;
    double numer = (1 + x2)*(6*x2 - 3) + 3 + x2*(2*x2 - 3);
    double denom = (1 + x2)*(6*x2 + 1) - 1 + x2*(2*x2 + 1);
    return sqrt(numer/(3*denom)) * C_LIGHT_CGS;
  }

  static void
  compute_pressure(body& particle){
    double rho = particle.getDensity();
    double Ye  = particle.getZbar()/particle.getAbar(); // TODO: shouldn't we be using electron fraction here?
    double P = pressure_from_rhoYe(rho, Ye);
    particle.setPressure(P);
  }

  /**
  * @brief      Compute sound speed for wd eos
  *
  * @param      particle
  */
  static void
  compute_soundspeed(body & particle) {
    double rho = particle.getDensity();
    double Ye  = particle.getZbar()/particle.getAbar(); // TODO: shouldn't we be using electron fraction here?
    double cs = soundspeed_from_rhoYe(rho, Ye);

#ifdef _DEBUG_EOS_
    if(cs != cs) {
      std::cout << "ERROR: speed of sound is NaN" << std::endl;
      std::cout << "Failed particle id: " << particle.id() << std::endl;
      std::cerr << "particle position: " << particle.coordinates() << std::endl;
      std::cerr << "particle velocity: " << particle.getVelocity() << std::endl;
      std::cerr << "particle acceleration: " << particle.getAcceleration()
                << std::endl;
      std::cerr << "smoothing length:  " << particle.radius() << std::endl;
      assert(false);
    }
#endif

    particle.setSoundspeed(cs);
  } // compute_soundspeed_wd

  /**
  * @brief      This is a zero-temperature equation of state, so temperature
  *             is decoupled from density and pressure; the function is thus
  *             empty
  *
  * @param      particle
  */
  static void
  compute_temperature(body & particle) {}

  /**
  * @brief      Compute specific internal energy
  *             TODO
  *
  * @param      particle
  */
  static void
  compute_internal_energy(body & particle) {
    // TODO: use piecewise polytropic approximation
    const double p = particle.getPressure(), rho = particle.getDensity();
    double u = 3./2. * p / rho;
    particle.setInternalenergy(u);
  }
}; // ...<eos_wd>

template<>
class eos_t<param::eos_ppt>{
  static double rho_thr;    // density threshold

public:
  static void init(body & particle){
    eos_t<param::eos_ppt>::rho_thr = param::ppt_density_thr;
    const double rho = particle.getDensity(),
                 P   = particle.getPressure();

    // In the piecewise-polytropic EOS, adiabatic invariant corresponds
    // to the first polytropic segment K1:
    //
    //  P(rho) = K1\rho^\Gamma1 + K2\rho^\Gamma2
    //
    double K1 = 0.0;
    if (rho < rho_thr) {
      K1 = P/pow(rho, poly_gamma);
    }
    else {
      double K2 = P/pow(rho, poly_gamma2);
      K1 = K2*pow(rho_thr, poly_gamma2 - poly_gamma);
    }
    particle.setAdiabatic(K1);
  }

  static void read_data(){}

  /**
  * @brief      Compute the pressure for piecewise-polytrope EOS
  * @param      particle
  */
  static void
  compute_pressure(body & particle) {
    const double rho = particle.getDensity(),
                 K1  = particle.getAdiabatic();
    double P = 0.0;
    if (rho < rho_thr) {
      P = K1*pow(rho, poly_gamma);
    }
    else {
      double K2 = K1*pow(rho_thr, poly_gamma - poly_gamma2);
      P = K2*pow(rho, poly_gamma2);
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
                 K1  = particle.getAdiabatic(),
                 gam = (rho < rho_thr ? poly_gamma : poly_gamma2);
    double soundspeed = 0.0;
    if (rho < rho_thr) {
      soundspeed = sqrt(K1*poly_gamma*pow(rho,poly_gamma - 2));
    }
    else {
      double K2 = K1*pow(rho_thr, poly_gamma - poly_gamma2);
      soundspeed = sqrt(K2*poly_gamma2*pow(rho,poly_gamma2 - 2));
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
                 K1  = particle.getAdiabatic();
    double eps = 0.0;
    if (rho < rho_thr) {
      eps = K1*pow(rho, poly_gamma - 1)/(poly_gamma - 1);
    }
    else {
      double K2 = K1*pow(rho_thr, poly_gamma - poly_gamma2);
      eps = K2*pow(rho,    poly_gamma2 - 1)/(poly_gamma2 - 1)
          - K2*pow(rho_thr,poly_gamma2 - 1)/(poly_gamma2 - 1)
          + K1*pow(rho_thr, poly_gamma - 1)/(poly_gamma - 1);
    }
    particle.setInternalenergy(eps);
  }

};

template<>
double eos_t<param::eos_ppt>::rho_thr;

template<>
class eos_t<param::eos_no_eos>{
public:
  static void init(body & particle){}
  static void read_data(){}
  static void compute_pressure(body& particle){}
  static void compute_soundspeed(body& particle){}
  static void compute_temperature(body& particle){}
  static void compute_internal_energy(body& particle){}
};

template<>
class eos_t<param::eos_wd_ideal_gas>{
public:
  static void init(body & particle){
    compute_internal_energy(particle);
  }
  static void read_data(){}

  static void compute_pressure(body& particle){
    const double
      rho = particle.getDensity(),
      eps = particle.getInternalenergy(),
      Ye  = particle.getZbar()/particle.getAbar();
    double P = (poly_gamma - 1.0)*rho*eps
             + eos_t<param::eos_wd>::pressure_from_rhoYe(rho, Ye);
    particle.setPressure(P);
  }

  /**
  * @brief      Compute sound speed for wd+ideal eos
  *
  * @param      particle
  */
  static void
  compute_soundspeed(body & particle) {
    const double
      rho = particle.getDensity(),
      eps = particle.getInternalenergy(),
      Ye  = particle.getZbar()/particle.getAbar();

    double cs2 = poly_gamma*(poly_gamma - 1)*eps
               + square(eos_t<param::eos_wd>::soundspeed_from_rhoYe(rho,Ye));
    double cs  = sqrt(cs2);

#ifdef _DEBUG_EOS_
    if(cs != cs) {
      std::cout << "ERROR: speed of sound is NaN" << std::endl;
      std::cout << "Failed particle id: " << particle.id() << std::endl;
      std::cerr << "particle position: " << particle.coordinates() << std::endl;
      std::cerr << "particle velocity: " << particle.getVelocity() << std::endl;
      std::cerr << "particle acceleration: " << particle.getAcceleration()
                << std::endl;
      std::cerr << "smoothing length:  " << particle.radius() << std::endl;
      assert(false);
    }
#endif

    particle.setSoundspeed(cs);
  }

  /**
  * @brief      Compute temperature, assuming ideal gas and
  *             fully ionized plasma
  *
  * @param      particle
  */
  static void
  compute_temperature(body & particle) {
    const double
        P = particle.getPressure(),
        rho = particle.getDensity(),
        abar = particle.getAbar(),
        zbar = particle.getZbar();
    double mu = abar*AMU/(zbar + 1.0);
    particle.setTemperature(P*mu/(rho*KBOL));
  }

  /**
  * @brief      Compute specific internal energy
  *             TODO
  *
  * @param      particle
  */
  static void
  compute_internal_energy(body & particle) {
    // TODO: check
    const double p = particle.getPressure(), rho = particle.getDensity();
    double u = 3./2. * p / rho;
    if (u < 0. ) log_one(error) << "u: " << u << std::endl;
    particle.setInternalenergy(u);
  }
};

// eos function types and pointers
typedef void (*compute_quantity_t)(body &);
typedef void (*read_data_t)();

#ifdef eos_type
#  define read_data            eos_t<eos_type>::read_data
#  define init                 eos_t<eos_type>::init
#  define compute_pressure     eos_t<eos_type>::compute_pressure
#  define compute_soundspeed   eos_t<eos_type>::compute_soundspeed
#  define compute_temperature  eos_t<eos_type>compute_temperature
#else
read_data_t read_data = nullptr;
compute_quantity_t init = nullptr;
compute_quantity_t compute_pressure = nullptr;
compute_quantity_t compute_soundspeed = nullptr;
compute_quantity_t compute_temperature = nullptr;
compute_quantity_t compute_internal_energy = nullptr;
#endif

/**
 * @brief  Installs the 'compute_pressure' and 'compute_soundspeed'
 *         function pointers, depending on the value of eos_type
 */
void
select() {
  using namespace param;

#ifndef eos_type
  switch(eos_type){
    case(eos_ideal):
      init = eos_t<eos_ideal>::init;
      read_data = eos_t<eos_ideal>::read_data;
      compute_pressure = eos_t<eos_ideal>::compute_pressure;
      compute_soundspeed = eos_t<eos_ideal>::compute_soundspeed;
      compute_temperature = eos_t<eos_ideal>::compute_temperature;
      compute_internal_energy = eos_t<eos_ideal>::compute_internal_energy;
      break;
    case(eos_polytropic):
      init = eos_t<eos_polytropic>::init;
      read_data = eos_t<eos_polytropic>::read_data;
      compute_pressure = eos_t<eos_polytropic>::compute_pressure;
      compute_soundspeed = eos_t<eos_polytropic>::compute_soundspeed;
      compute_temperature = eos_t<eos_polytropic>::compute_temperature;
      compute_internal_energy = eos_t<eos_polytropic>::compute_internal_energy;
      break;
    case(eos_wd):
      init = eos_t<eos_wd>::init;
      read_data = eos_t<eos_wd>::read_data;
      compute_pressure = eos_t<eos_wd>::compute_pressure;
      compute_soundspeed = eos_t<eos_wd>::compute_soundspeed;
      compute_temperature = eos_t<eos_wd>::compute_temperature;
      compute_internal_energy = eos_t<eos_wd>::compute_internal_energy;
      break;
    case(eos_ppt):
      init = eos_t<eos_ppt>::init;
      read_data = eos_t<eos_ppt>::read_data;
      compute_pressure = eos_t<eos_ppt>::compute_pressure;
      compute_soundspeed = eos_t<eos_ppt>::compute_soundspeed;
      compute_temperature = eos_t<eos_ppt>::compute_temperature;
      compute_internal_energy = eos_t<eos_ppt>::compute_internal_energy;
      break;
    case(eos_no_eos):
      init = eos_t<eos_no_eos>::init;
      read_data = eos_t<eos_no_eos>::read_data;
      compute_pressure = eos_t<eos_no_eos>::compute_pressure;
      compute_soundspeed = eos_t<eos_no_eos>::compute_soundspeed;
      compute_temperature = eos_t<eos_no_eos>::compute_temperature;
      compute_internal_energy = eos_t<eos_no_eos>::compute_internal_energy;
      break;
    case(eos_stellar_collapse):
      init = eos_t<eos_stellar_collapse>::init;
      read_data = eos_t<eos_stellar_collapse>::read_data;
      compute_pressure = eos_t<eos_stellar_collapse>::compute_pressure;
      compute_soundspeed = eos_t<eos_stellar_collapse>::compute_soundspeed;
      compute_temperature = eos_t<eos_stellar_collapse>::compute_temperature;
      compute_internal_energy = eos_t<eos_stellar_collapse>::compute_internal_energy;
      break;
    case(eos_wd_ideal_gas):
      init = eos_t<eos_wd_ideal_gas>::init;
      read_data = eos_t<eos_wd_ideal_gas>::read_data;
      compute_pressure = eos_t<eos_wd_ideal_gas>::compute_pressure;
      compute_soundspeed = eos_t<eos_wd_ideal_gas>::compute_soundspeed;
      compute_temperature = eos_t<eos_wd_ideal_gas>::compute_temperature;
      compute_internal_energy = eos_t<eos_wd_ideal_gas>::compute_internal_energy;
      break;
    case(eos_helmholtz):
      init = eos_t<eos_helmholtz>::init;
      read_data = eos_t<eos_helmholtz>::read_data;
      compute_pressure = eos_t<eos_helmholtz>::compute_pressure;
      compute_soundspeed = eos_t<eos_helmholtz>::compute_soundspeed;
      compute_temperature = eos_t<eos_helmholtz>::compute_temperature;
      compute_internal_energy = eos_t<eos_helmholtz>::compute_internal_energy;
      break;
    default:
      init = nullptr;
      read_data = nullptr;
      compute_pressure = nullptr;
      compute_soundspeed = nullptr;
      compute_temperature = nullptr;
      std::cerr<<"Undefined eos type"<<std::endl;
      MPI_Finalize();
      exit(0);
      break;
  }
#endif // eos_type
} // select

/**
 * @brief      set uniform average atomic weight (abar) and proton number (zbar)
 *             for the particle, from initial abar and zbar params
 *             TODO: read from species file
 *
 * @param      particle
 *
 * @uses       initial_abar     global parameter
 * @uses       initial_zbar     global parameter
 */
void
initialize_abarzbar(
  body& particle)
{
  using namespace param;
  particle.setAbar(initial_abar);
  particle.setZbar(initial_zbar);
} // initialize_abarzbar

} // namespace eos

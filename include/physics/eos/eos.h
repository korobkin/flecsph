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

constexpr double square(const double& x){
  return ((x) * (x));
}
constexpr double cube(const double& x){
  return ((x) * (x) * (x));
}
constexpr double quartic(const double& x){
  return ((x) * (x) * (x) * (x));
}

using namespace param;

template<>
class eos_t<param::eos_ideal>{
public:
  static void init(body & particle){}

  static void read_data(){}

  static void compute_pressure(body& particle){
    using namespace param;
    double pressure =
      (poly_gamma - 1.0) * particle.getDensity() * particle.getInternalenergy();
    particle.setPressure(pressure);
  }
  /**
  * @brief      Compute sound speed for ideal fluid or polytropic eos
  *
  * @param      particle
  */
  static void compute_soundspeed(body & particle) {
    using namespace param;
    double soundspeed =
      sqrt(poly_gamma * particle.getPressure() / particle.getDensity());
    particle.setSoundspeed(soundspeed);
  }

  /**
  * @brief      Compute temperature via ideal gas formula
  *
  * @param      particle
  */
  static void compute_temperature(body & particle) {
    const double abar = particle.getAbar(),
                 P = particle.getPressure(),
                 rho = particle.getDensity();
    double T = AMU*abar*P/(rho*KBOL);
    particle.setTemperature(T);
  }
}; // ...<eos_ideal>

template<>
class eos_t<param::eos_polytropic>{
public:
  /**
  * @brief      Equation-of-state intializer:
  *             computes missing quantities etc.
  * @param      particle
  */
  static void init(body & particle){
    using namespace param;
    double K = particle.getPressure() / pow(particle.getDensity(), poly_gamma);
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
    using namespace param;
    double pressure =
      particle.getAdiabatic() * pow(particle.getDensity(), poly_gamma);
    particle.setPressure(pressure);
  }

  /**
  * @brief      Compute sound speed for ideal fluid or polytropic eos
  *
  * @param      particle
  */
  static void compute_soundspeed(body & particle) {
    using namespace param;
    double soundspeed =
      sqrt(poly_gamma * particle.getPressure() / particle.getDensity());
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

  static constexpr double A_wd = 6.00288e22;
  static constexpr double B_wd_nm = 9.81011e5;

public:
  static void init(body & particle){
    setInternalenergy(particle);
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

private:
  static void setInternalenergy(body & particle){
    const double p = particle.getPressure(), rho = particle.getDensity();
    double u = 3./2. * p / rho;
    if (u < 0. ) log_one(error) << "u: " << u << std::endl;
    particle.setInternalenergy(u);
  }
}; // ...<eos_wd>

template<>
class eos_t<param::eos_ppt>{
public:
  static void init(body & particle){}

  static void read_data(){}

  /**
  * @brief      Compute the pressure for piecewise-polytrope EOS
  * @param      particle
  */
  static void
  compute_pressure(body & particle) {
    using namespace param;
    // TODO : transient density might be parametrized
    //       or determining automatically.
    //       But here I put certian value that is
    //       reasonalbe transition between
    //       relativisitc and non-relativistic regimes
    const double transition_density = 5e+14;

    if(particle.getDensity() <= transition_density) {
      double pressure =
        particle.getAdiabatic() * pow(particle.getDensity(), poly_gamma);
      particle.setPressure(pressure);
    }
    else {
      double pressure =
        (particle.getAdiabatic() * pow(transition_density, poly_gamma) /
          pow(transition_density, poly_gamma2)) *
        pow(particle.getDensity(), poly_gamma2);
      particle.setPressure(pressure);
    }
  } // compute_pressure_ppt

  /**
  * @brief      Compute sound speed for piecewise polytropic eos
  *
  * @param      particle
  */
  static void
  compute_soundspeed(body & particle) {
    using namespace param;
    double density_transition = 500000000000000;
    double density = particle.getDensity();
    if(density <= density_transition) {
      double soundspeed =
        sqrt(poly_gamma * particle.getPressure() / density);
      particle.setSoundspeed(soundspeed);
    }
    else {
      double soundspeed =
        sqrt(poly_gamma2 * particle.getPressure() / density);
      particle.setSoundspeed(soundspeed);
    }
  } // compute_soundspeed_ppt

  /**
  * @brief      Empty function because EOS is temperature-independent
  *
  * @param      particle
  */
  static void
  compute_temperature(body&){}
};

template<>
class eos_t<param::eos_no_eos>{
public:
  static void init(body & particle){}
  static void read_data(){}
  static void compute_pressure(body& particle){}
  static void compute_soundspeed(body& particle){}
  static void compute_temperature(body& particle){}
};

template<>
class eos_t<param::eos_wd_ideal_gas>{
public:
  static void init(body & particle){
    setInternalenergy(particle);
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

    double P_ideal = (poly_gamma - 1.0)*rho*eps;
    double cs2 = poly_gamma*P_ideal/rho 
               + square(soundspeed_from_rhoYe(rho,Ye);

    particle.setSoundspeed(sqrt(cs2));
  }

  /**
  * @brief      Compute temperature, assuming ideal gas and 
  *             fully ionized plasma
  *
  * @param      particle
  */
  static void
  compute_temperature(body & particle) {
    const double p = particle.getPressure(), rho = particle.getDensity(),
              abar = particle.getAbar(), zbar = particle.getZbar();
    double mu = abar * AMU / (zbar + 1.0);
    double T = mu * p / (rho * KBOL);
    particle.setTemperature(T);
  }
  private:
  static void setInternalenergy(body & particle){
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
      break;
    case(eos_polytropic):
      init = eos_t<eos_polytropic>::init;
      read_data = eos_t<eos_polytropic>::read_data;
      compute_pressure = eos_t<eos_polytropic>::compute_pressure;
      compute_soundspeed = eos_t<eos_polytropic>::compute_soundspeed;
      compute_temperature = eos_t<eos_polytropic>::compute_temperature;
      break;
    case(eos_wd):
      init = eos_t<eos_wd>::init;
      read_data = eos_t<eos_wd>::read_data;
      compute_pressure = eos_t<eos_wd>::compute_pressure;
      compute_soundspeed = eos_t<eos_wd>::compute_soundspeed;
      compute_temperature = eos_t<eos_wd>::compute_temperature;
      break;
    case(eos_ppt):
      init = eos_t<eos_ppt>::init;
      read_data = eos_t<eos_ppt>::read_data;
      compute_pressure = eos_t<eos_ppt>::compute_pressure;
      compute_soundspeed = eos_t<eos_ppt>::compute_soundspeed;
      compute_temperature = eos_t<eos_ppt>::compute_temperature;
      break;
    case(eos_no_eos):
      init = eos_t<eos_no_eos>::init;
      read_data = eos_t<eos_no_eos>::read_data;
      compute_pressure = eos_t<eos_no_eos>::compute_pressure;
      compute_soundspeed = eos_t<eos_no_eos>::compute_soundspeed;
      compute_temperature = eos_t<eos_no_eos>::compute_temperature;
      break;
    case(eos_stellar_collapse):
      init = eos_t<eos_stellar_collapse>::init;
      read_data = eos_t<eos_stellar_collapse>::read_data;
      compute_pressure = eos_t<eos_stellar_collapse>::compute_pressure;
      compute_soundspeed = eos_t<eos_stellar_collapse>::compute_soundspeed;
      compute_temperature = eos_t<eos_stellar_collapse>::compute_temperature;
      break;
    case(eos_wd_ideal_gas):
      init = eos_t<eos_wd_ideal_gas>::init;
      read_data = eos_t<eos_wd_ideal_gas>::read_data;
      compute_pressure = eos_t<eos_wd_ideal_gas>::compute_pressure;
      compute_soundspeed = eos_t<eos_wd_ideal_gas>::compute_soundspeed;
      compute_temperature = eos_t<eos_wd_ideal_gas>::compute_temperature;
      break;
    case(eos_helmholtz):
      init = eos_t<eos_helmholtz>::init;
      read_data = eos_t<eos_helmholtz>::read_data;
      compute_pressure = eos_t<eos_helmholtz>::compute_pressure;
      compute_soundspeed = eos_t<eos_helmholtz>::compute_soundspeed;
      compute_temperature = eos_t<eos_helmholtz>::compute_temperature;
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

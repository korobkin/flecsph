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
  static void init(){}

  static void read_data(){}

  static void compute_pressure(body& source){
    using namespace param;
    double pressure =
      (poly_gamma - 1.0) * source.getDensity() * source.getInternalenergy();
    source.setPressure(pressure);
  }
  /**
  * @brief      Compute sound speed for ideal fluid or polytropic eos
  * From CES-Seminar 13/14 - Smoothed Particle Hydrodynamics
  *
  * @param      srch  The source's body holder
  */
  static void compute_soundspeed(body & source) {
    using namespace param;
    double soundspeed =
      sqrt(poly_gamma * source.getPressure() / source.getDensity());
    source.setSoundspeed(soundspeed);
  }

  /**
  * @brief      Compute temperature via ideal gas WD
  *             TODO: double-check formula [???]
  *
  * @param      srch  The source's body holder
  */
  static void compute_temperature(body & source) {
    const double abar = source.getAbar(), zbar = source.getZbar(),
                 P = source.getPressure(), rho = source.getDensity(),
                 Ye = source.getElectronfraction();
    double mu = abar * (AMU + Ye * ME) / (zbar + 1.0); // ???
    double T = mu * P / (rho * KBOL);
    source.setTemperature(T);
  } // compute_temperature_ideal
};

template<>
class eos_t<param::eos_polytropic>{
public:
  **
  * @brief      Equation-of-state intializer:
  *             computes missing quantities etc.
  * @param      srch  The source's body holder
  *
  static void init(){
    using namespace param;
    double K = source.getPressure() / pow(source.getDensity(), poly_gamma);
    source.setAdiabatic(K);
    return;
  }

  static void read_data(){}

  **
  * @brief      Compute pressure from density using polytrope
  *             P(\rho) = A*\rho^\Gamma
  *
  * @param      srch  The source's body holder
  *
  static void compute_pressure(body & source) {
    using namespace param;
    double pressure =
      source.getAdiabatic() * pow(source.getDensity(), poly_gamma);
    source.setPressure(pressure);
  }

  **
  * @brief      Compute sound speed for ideal fluid or polytropic eos
  * From CES-Seminar 13/14 - Smoothed Particle Hydrodynamics
  *
  * @param      srch  The source's body holder
  *
  static void compute_soundspeed(body & source) {
    using namespace param;
    double soundspeed =
      sqrt(poly_gamma * source.getPressure() / source.getDensity());
    source.setSoundspeed(soundspeed);
  }

  static void
  compute_temperature(body& source){
    assert(false&&"Not implemented");
  }
};

template<>
class eos_t<param::eos_wd>{
public:
  static void init(){}

  static void read_data(){}

  static void compute_pressure(body& source){
    double Ye   = source.getZbar()/source.getAbar();
    double A_wd = 6.00288e22;
    double B_wd = 9.81011e5 / Ye;

    double x_wd = pow(source.getDensity() / B_wd, 1.0 / 3.0);
    double pressure =
      A_wd *
      (x_wd * (2.0 * square(x_wd) - 3.0) * sqrt(square(x_wd) + 1.0) + 3.0 * asinh(x_wd));
    source.setPressure(pressure);
  }

  /**
  * @brief      Compute sound speed for wd eos
  *
  * @param      srch  The source's body holder
  */
  static void
  compute_soundspeed(body & source) {
    using namespace param;
    double Ye   = source.getZbar()/source.getAbar();
    double A_wd = 6.00288e22;
    double B_wd = 9.81011e5 / Ye;
    double x_wd = pow(source.getDensity() / B_wd, 1. / 3.);

    double sterm = sqrt(1.0 + square(x_wd));
    double numer = 3.0 / sterm + sterm * (6.0 * square(x_wd) - 3.0) +
                  square(x_wd) * (2.0 * square(x_wd) - 3.0) / sterm;
    double denom = -1.0 / sterm + sterm * (1.0 + 6.0 * square(x_wd)) +
                  square(x_wd) * (1.0 + 2.0 * square(x_wd)) / sterm;

    double soundspeed = sqrt(numer / (3.0 * denom)) * C_LIGHT_CGS;

    if(not(numer / denom > 0)) {
      std::cout << "speed of sounds is not a real number: "
                << "numer/denom = " << numer / denom << std::endl;
      std::cout << "Failed particle id: " << source.id() << std::endl;
      std::cerr << "particle position: " << source.coordinates() << std::endl;
      std::cerr << "particle velocity: " << source.getVelocity() << std::endl;
      std::cerr << "particle acceleration: " << source.getAcceleration()
                << std::endl;
      std::cerr << "smoothing length:  " << source.radius() << std::endl;
      assert(false);
    }

    // double numer = 8.*source.getDensity()*x_wd - 3.*B_wd;
    // double deno = 3*B_wd*B_wd*x_wd*x_wd*sqrt(x_wd*x_wd+1);

    // double soundspeed = A_wd*(numer/deno +
    //                          x_wd/(3.*source.getDensity()
    //                                *sqrt(1-x_wd*x_wd)));
    source.setSoundspeed(soundspeed);
  } // compute_soundspeed_wd

  /**
  * @brief      Compute temperature via ideal gas in C/O WD
  *             TODO: parameterize abar, zbar; double-check formula [???]
  *
  * @param      srch  The source's body holder
  */
  static void
  compute_temperature(body & source) {
    const double abar = source.getAbar(), zbar = source.getZbar(),
                 P = source.getPressure(), rho = source.getDensity(),
                 Ye = zbar/abar;
    double mu = abar * (AMU + Ye * ME) / (zbar + 1.0); // ???
    double T = mu * P / (rho * KBOL);
    source.setTemperature(T);
  } // compute_temperature_ideal

};

template<>
class eos_t<param::eos_ppt>{
public:
  static void init(){}

  static void read_data(){}

  /**
  * @brief      Compute the pressure for piecewise-polytrope EOS
  * @param      srch  The source's body holder
  */
  static void
  compute_pressure(body & source) {
    using namespace param;
    // TODO : transient density might be parametrized
    //       or determining automatically.
    //       But here I put certian value that is
    //       reasonalbe transition between
    //       relativisitc and non-relativistic regimes
    const double transition_density = 5e+14;

    if(source.getDensity() <= transition_density) {
      double pressure =
        source.getAdiabatic() * pow(source.getDensity(), poly_gamma);
      source.setPressure(pressure);
    }
    else {
      double pressure =
        (source.getAdiabatic() * pow(transition_density, poly_gamma) /
          pow(transition_density, poly_gamma2)) *
        pow(source.getDensity(), poly_gamma2);
      source.setPressure(pressure);
    }
  } // compute_pressure_ppt

  /**
  * @brief      Compute sound speed for piecewise polytropic eos
  *
  * @param      srch  The source's body holder
  */
  static void
  compute_soundspeed(body & source) {
    using namespace param;
    double density_transition = 500000000000000;
    double density = source.getDensity();
    if(density <= density_transition) {
      double soundspeed =
        sqrt(poly_gamma * source.getPressure() / density);
      source.setSoundspeed(soundspeed);
    }
    else {
      double soundspeed =
        sqrt(poly_gamma2 * source.getPressure() / density);
      source.setSoundspeed(soundspeed);
    }
  } // compute_soundspeed_ppt

  static void
  compute_temperature(body&){}
};

template<>
class eos_t<param::eos_no_eos>{
public:
  static void init(){}
  static void read_data(){}
  static void compute_pressure(body& source){}
  static void compute_soundspeed(body& source){}
  static void compute_temperature(body& source){}
};

template<>
class eos_t<param::eos_pure_gravitation>{
public:
  static void init(){}
  static void read_data(){}
  static void compute_pressure(body& source){}
  /**
  * @brief      Compute sound speed for ideal fluid or polytropic eos
  * From CES-Seminar 13/14 - Smoothed Particle Hydrodynamics
  *
  * @param      srch  The source's body holder
  */
  static void
  compute_soundspeed(body & source) {
    using namespace param;
    double soundspeed =
      sqrt(poly_gamma * source.getPressure() / source.getDensity());
    source.setSoundspeed(soundspeed);
  }
  static void compute_temperature(body& source){}
};

template<>
class eos_t<param::eos_wd_ideal_gas>{
public:
  static void init(){}
  static void read_data(){}
  
  static void compute_pressure(body& source){
    double density = source.getDensity();
    double pressure =
      (poly_gamma - 1.0) * density * source.getInternalenergy();
    double Ye   = source.getZbar()/source.getAbar();
    double A_wd = 6.00288e22;
    double B_wd = 9.81011e5 / Ye;

    double x_wd = pow(density / B_wd, 1.0 / 3.0);
    pressure +=
      A_wd *
      (x_wd * (2.0 * square(x_wd) - 3.0) * sqrt(square(x_wd) + 1.0) + 3.0 * asinh(x_wd));
    source.setPressure(pressure);
  }

  /**
  * @brief      Compute sound speed for wd eos
  *
  * @param      srch  The source's body holder
  */
  static void
  compute_soundspeed(body & source) {
    double density = source.getDensity();
    double ideal_pressure =
      (poly_gamma - 1.0) * density * source.getInternalenergy();
    double soundspeed = poly_gamma * ideal_pressure / density;

    double Ye   = source.getZbar()/source.getAbar();
    double A_wd = 6.00288e22;
    double B_wd = 9.81011e5 / Ye;
    double x_wd = pow(density / B_wd, 1. / 3.);

    double sterm = sqrt(1.0 + square(x_wd));
    double numer = 3.0 / sterm + sterm * (6.0 * square(x_wd) - 3.0) +
                  square(x_wd) * (2.0 * square(x_wd) - 3.0) / sterm;
    double denom = -1.0 / sterm + sterm * (1.0 + 6.0 * square(x_wd)) +
                  square(x_wd) * (1.0 + 2.0 * square(x_wd)) / sterm;

    soundspeed += (numer / (3.0 * denom)) * square(C_LIGHT_CGS);

    if(not(numer / denom > 0)) {
      std::cout << "speed of sounds is not a real number: "
                << "numer/denom = " << numer / denom << std::endl;
      std::cout << "Failed particle id: " << source.id() << std::endl;
      std::cerr << "particle position: " << source.coordinates() << std::endl;
      std::cerr << "particle velocity: " << source.getVelocity() << std::endl;
      std::cerr << "particle acceleration: " << source.getAcceleration()
                << std::endl;
      std::cerr << "smoothing length:  " << source.radius() << std::endl;
      assert(false);
    }

    // double numer = 8.*source.getDensity()*x_wd - 3.*B_wd;
    // double deno = 3*B_wd*B_wd*x_wd*x_wd*sqrt(x_wd*x_wd+1);

    // double soundspeed = A_wd*(numer/deno +
    //                          x_wd/(3.*source.getDensity()
    //                                *sqrt(1-x_wd*x_wd)));
    source.setSoundspeed(sqrt(soundspeed));
  } // compute_soundspeed_wd_ideal

  /**
  * @brief      Compute temperature via ideal gas in C/O WD
  *             TODO: double-check formula [???]
  *
  * @param      srch  The source's body holder
  */
  static void
  compute_temperature(body & source) {
    const double P = source.getPressure(), rho = source.getDensity(),
              abar = source.getAbar(), zbar = source.getZbar();
    double Ye = zbar/abar;
    double mu = abar * (AMU + Ye * ME) / (zbar + 1.0); // ???
    double T = mu * P / (rho * KBOL);
    source.setTemperature(T);
  } // compute_temperature_ideal

};

// eos function types and pointers
typedef void (*compute_quantity_t)(body &);
typedef void (*data_read)();

#ifdef eos_type
constexpr data_read read_data = eos_t<eos_type>::read_data;
constexpr compute_quantity_t init = eos_t<eos_type>::init;
constexpr compute_quantity_t compute_pressure = eos_t<eos_type>::compute_pressure;
constexpr compute_quantity_t compute_soundspeed = eos_t<eos_type>::compute_soundspeed;
constexpr compute_quantity_t compute_temperature = eos_t<eos_type>compute_temperature;
#else
data_read read_data = nullptr;
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
    case(eos_pure_gravitation):
      init = eos_t<eos_pure_gravitation>::init;
      read_data = eos_t<eos_pure_gravitation>::read_data;
      compute_pressure = eos_t<eos_pure_gravitation>::compute_pressure;
      compute_soundspeed = eos_t<eos_pure_gravitation>::compute_soundspeed;
      compute_temperature = eos_t<eos_pure_gravitation>::compute_temperature;
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
 * @brief      set the abar and zbar for each particle
 *             TODO: read from species file
 *
 * @param      src       The source particle
 *
 * @return
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

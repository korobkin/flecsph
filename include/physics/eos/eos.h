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
#include "eos_stellar_collapse.h"

namespace eos {

constexpr double square(const double& x){
  return ((x) * (x));
}
constexpr double cube(const double& x){
  return ((x) * (x) * (x));
}

using namespace param;

template<>
class eos_t<param::eos_ideal>{
public: 
  static void init(body&){}
  
  static void compute_pressure(body & particle){
    using namespace param;
    double pressure =
      (poly_gamma - 1)*particle.getDensity()*particle.getInternalenergy();
    particle.setPressure(pressure);
  }
  /**
  * @brief      Compute sound speed for ideal fluid or polytropic eos
  * From CES-Seminar 13/14 - Smoothed Particle Hydrodynamics
  *
  * @param      particle  The particle in question
  */
  static void compute_soundspeed(body & particle) {
    using namespace param;
    double soundspeed =
      sqrt(poly_gamma * particle.getPressure() / particle.getDensity());
    particle.setSoundspeed(soundspeed);
  }

  /**
  * @brief      Compute temperature via ideal gas in C/O WD
  *             TODO: parameterize abar, zbar; double-check formula [???]
  *
  * @param      particle  The particle in question
  */
  static void compute_temperature(body & particle) {
    const double kB = 1.3806505e-16, // [erg/K]
      abar = 12.0, // [mol/g] molar mass of Carbon-12
      zbar = 6.0, // proton number for C
      amu = 1.66053906660e-24, // [g] a.m.u.
      me = 9.10938356e-28; // [g] electron mass
    const double P = particle.getPressure(), rho = particle.getDensity(),
                Ye = particle.getElectronfraction();
    double mu = abar * (amu + Ye * me) / (zbar + 1.0); // ???
    double T = mu * P / (rho * kB);
    particle.setTemperature(T);
  } // compute_temperature_ideal
};

template<>
class eos_t<param::eos_polytropic>{
public: 
    /**
  * @brief      Equation-of-state intializer:
  *             computes missing quantities etc.
  * @param      particle  The particle in question
  */
  static void init(body & particle){
    using namespace param;
    double K = particle.getPressure()/pow(particle.getDensity(),poly_gamma);
    particle.setAdiabatic(K);
    return;
  }
  /**
  * @brief      Compute pressure from density using polytrope
  *             P(\rho) = A*\rho^\Gamma
  *
  * @param      particle  The particle in question
  */
  static void compute_pressure(body & particle) {
    using namespace param;
    double pressure =
      particle.getAdiabatic()*pow(particle.getDensity(),poly_gamma);
    particle.setPressure(pressure);
  }

  /**
  * @brief      Compute sound speed for ideal fluid or polytropic eos
  * From CES-Seminar 13/14 - Smoothed Particle Hydrodynamics
  *
  * @param      particle  The particle in question
  */
  static void compute_soundspeed(body & particle) {
    using namespace param;
    double soundspeed =
      sqrt(poly_gamma * particle.getPressure() / particle.getDensity());
    particle.setSoundspeed(soundspeed);
  }

  static void 
  compute_temperature(body & particle){
    assert(false&&"Not implemented");
  }
};


template<>
class eos_t<param::eos_wd>{
public: 
  static void init(body&){}
  static void compute_pressure(body & particle){
    double Ye = 0.5;
    double A_wd = 6.00288e22;
    double B_wd = 9.81011e5 / Ye;

    double x_wd = pow((particle.getDensity())/B_wd, 1./3.);
    double pressure =
      A_wd*(x_wd*(2*square(x_wd) - 3)*sqrt(square(x_wd) + 1) + 3*asinh(x_wd));
    particle.setPressure(pressure);
  }

  /**
  * @brief      Compute sound speed for wd eos
  *
  * @param      particle  The particle in question
  */
  static void
  compute_soundspeed(body & particle) {
    using namespace param;
    double Ye = 0.5;
    double A_wd = 6.00288e22;
    double B_wd = 9.81011e5 / Ye;
    double x_wd = pow((particle.getDensity()) / B_wd, 1. / 3.);
    double cc = 29979245800.0;

    double sterm = sqrt(1.0 + square(x_wd));
    double numer = 3.0 / sterm + sterm * (6.0 * square(x_wd) - 3.0) +
                  square(x_wd) * (2.0 * square(x_wd) - 3.0) / sterm;
    double denom = -1.0 / sterm + sterm * (1.0 + 6.0 * square(x_wd)) +
                  square(x_wd) * (1.0 + 2.0 * square(x_wd)) / sterm;

    double soundspeed = sqrt(numer / (3.0 * denom)) * cc;

    if(not(numer / denom > 0)) {
      std::cout << "speed of sounds is not a real number: "
                << "numer/denom = " << numer / denom << std::endl;
      std::cout << "Failed particle id: " << particle.id() << std::endl;
      std::cerr << "particle position: " << particle.coordinates() << std::endl;
      std::cerr << "particle velocity: " << particle.getVelocity() << std::endl;
      std::cerr << "particle acceleration: " << particle.getAcceleration()
                << std::endl;
      std::cerr << "smoothing length:  " << particle.radius() << std::endl;
      assert(false);
    }

    // double numer = 8.*particle.getDensity()*x_wd - 3.*B_wd;
    // double deno = 3*B_wd*B_wd*x_wd*x_wd*sqrt(x_wd*x_wd+1);

    // double soundspeed = A_wd*(numer/deno +
    //                          x_wd/(3.*particle.getDensity()
    //                                *sqrt(1-x_wd*x_wd)));
    particle.setSoundspeed(soundspeed);
  } // compute_soundspeed_wd

  /**
  * @brief      Compute temperature via ideal gas in C/O WD
  *             TODO: parameterize abar, zbar; double-check formula [???]
  *
  * @param      particle  The particle in question
  */
  static void
  compute_temperature(body & particle) {
    const double kB = 1.3806505e-16, // [erg/K]
      abar = 12.0, // [mol/g] molar mass of Carbon-12
      zbar = 6.0, // proton number for C
      amu = 1.66053906660e-24, // [g] a.m.u.
      me = 9.10938356e-28; // [g] electron mass
    const double P = particle.getPressure(), rho = particle.getDensity(),
                Ye = particle.getElectronfraction();
    double mu = abar * (amu + Ye * me) / (zbar + 1.0); // ???
    double T = mu * P / (rho * kB);
    particle.setTemperature(T);
  } // compute_temperature_ideal

};

template<>
class eos_t<param::eos_ppt>{
public: 
  static void init(body&){}
  /**
  * @brief      Compute the pressure for piecewise-polytrope EOS
  * @param      particle  The particle in question
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
  * @param      particle  The particle in question
  */
  static void
  compute_soundspeed(body & particle) {
    using namespace param;
    double density_transition = 500000000000000;
    if(particle.getDensity() <= density_transition) {
      double soundspeed =
        sqrt(poly_gamma * particle.getPressure() / particle.getDensity());
      particle.setSoundspeed(soundspeed);
    }
    else {
      double soundspeed =
        sqrt(poly_gamma2 * particle.getPressure() / particle.getDensity());
      particle.setSoundspeed(soundspeed);
    }
  } // compute_soundspeed_ppt
  static void 
  compute_temperature(body&){}
}; 

template<>
class eos_t<param::eos_no_eos>{
public: 
  static void init(body&){}
  static void compute_pressure(body&){}
  static void compute_soundspeed(body&){}
  static void compute_temperature(body&){}
};

template<>
class eos_t<param::eos_pure_gravitation>{
public: 
  static void init(body&){}
  static void compute_pressure(body&){}
  /**
  * @brief      Compute sound speed for ideal fluid or polytropic eos
  * From CES-Seminar 13/14 - Smoothed Particle Hydrodynamics
  *
  * @param      particle  The particle in question
  */
  static void
  compute_soundspeed(body & particle) {
    using namespace param;
    double soundspeed =
      sqrt(poly_gamma * particle.getPressure() / particle.getDensity());
    particle.setSoundspeed(soundspeed);
  }
  static void compute_temperature(body& particle){}
};

// eos function types and pointers
typedef void (*compute_quantity_t)(body &);

#ifdef eos_type 
# define eos_init             eos_t<param::eos_type>::init
# define compute_pressure     eos_t<param::eos_type>::compute_pressure
# define compute_soundspeed   eos_t<param::eos_type>::compute_soundspeed
# define compute_temperature  eos_t<param::eos_type>::compute_temperature
#else 
compute_quantity_t eos_init = nullptr; 
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
      eos_init = eos_t<eos_ideal>::init; 
      compute_pressure = eos_t<eos_ideal>::compute_pressure;
      compute_soundspeed = eos_t<eos_ideal>::compute_soundspeed;
      compute_temperature = eos_t<eos_ideal>::compute_temperature;
      break;
    case(eos_polytropic): 
      eos_init = eos_t<eos_polytropic>::init; 
      compute_pressure = eos_t<eos_polytropic>::compute_pressure;
      compute_soundspeed = eos_t<eos_polytropic>::compute_soundspeed;
      compute_temperature = eos_t<eos_polytropic>::compute_temperature;
      break;
    case(eos_wd): 
      eos_init = eos_t<eos_wd>::init; 
      compute_pressure = eos_t<eos_wd>::compute_pressure;
      compute_soundspeed = eos_t<eos_wd>::compute_soundspeed;
      compute_temperature = eos_t<eos_wd>::compute_temperature;
      break;
    case(eos_ppt): 
      eos_init = eos_t<eos_ppt>::init; 
      compute_pressure = eos_t<eos_ppt>::compute_pressure;
      compute_soundspeed = eos_t<eos_ppt>::compute_soundspeed;
      compute_temperature = eos_t<eos_ppt>::compute_temperature;
      break;
    case(eos_no_eos): 
      eos_init = eos_t<eos_no_eos>::init; 
      compute_pressure = eos_t<eos_no_eos>::compute_pressure;
      compute_soundspeed = eos_t<eos_no_eos>::compute_soundspeed;
      compute_temperature = eos_t<eos_no_eos>::compute_temperature;
      break;
    case(eos_pure_gravitation): 
      eos_init = eos_t<eos_pure_gravitation>::init; 
      compute_pressure = eos_t<eos_pure_gravitation>::compute_pressure;
      compute_soundspeed = eos_t<eos_pure_gravitation>::compute_soundspeed;
      compute_temperature = eos_t<eos_pure_gravitation>::compute_temperature;
      break;
    case(eos_stellar_collapse): 
      eos_init = eos_t<eos_stellar_collapse>::init; 
      compute_pressure = eos_t<eos_stellar_collapse>::compute_pressure;
      compute_soundspeed = eos_t<eos_stellar_collapse>::compute_soundspeed;
      compute_temperature = eos_t<eos_stellar_collapse>::compute_temperature;
      break;
    default: 
      eos_init = nullptr; 
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

} // namespace eos

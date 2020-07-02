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
class eos_t<param::eos_polytropic>{

public:
  /**
  * @brief      Compute adiabatic invariant from density and pressure
  *
  * @param      rho   density
  * @param      P     pressure
  */
  static inline double
  adiabatic_given_rhoP (const double rho, const double P) {
    return P/pow(rho,poly_gamma);
  }

  /**
  * @brief      Compute adiabatic invariant from density and pressure
  *
  * @param      particle
  */
  static void
  compute_adiabatic(body & particle){
    const double rho = particle.getDensity(),
                 P = particle.getPressure();
    double K = adiabatic_given_rhoP(rho, P);
    particle.setAdiabatic(K);
  }

  /**
  * @brief      Initialized adiabatic invariant from initial conditions
  *
  * @param      particle
  */
  static void init(body & particle){
    compute_adiabatic(particle);
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
    double soundspeed = sqrt(K*poly_gamma*pow(rho, poly_gamma - 2.));
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
    double eps = K*pow(rho, poly_gamma - 1.)/(poly_gamma - 1.);
    particle.setInternalenergy(eps);
  }

}; // ...<eos_polytropic>


template<>
class eos_t<param::eos_ideal>{

public:
  /**
  * @brief      Initialize missing thermodynamic quantities
  *
  * @param      particle
  */
  static void init(body & particle){
    eos_t<param::eos_polytropic>::compute_adiabatic(particle);
  }

  static void read_data(){}

  /**
  * @brief      Computes pressure using the density and internal energy
  *
  * @param      particle
  */
  static void
  compute_pressure(body& particle){
    double pressure =
      (poly_gamma - 1.) * particle.getDensity() * particle.getInternalenergy();
    particle.setPressure(pressure);
  }

  /**
  * @brief      Sound speed from specific internal energy
  *
  * @param      particle
  */
  static void
  compute_soundspeed(body & particle) {
    const double eps = particle.getInternalenergy();
    double soundspeed = sqrt(poly_gamma*(poly_gamma - 1.)*eps);
    particle.setSoundspeed(soundspeed);
  }

  /**
  * @brief      Compute temperature via ideal gas formula
  *             Uses abar and specific internal energy
  *
  * @param      particle
  */
  static void
  compute_temperature(body & particle) {
    const double abar = particle.getAbar(),
                 eps  = particle.getInternalenergy();
    double T = AMU/KBOL*abar*(poly_gamma - 1.)*eps;
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
    double eps = K*pow(rho, poly_gamma - 1.)/(poly_gamma - 1.);
    particle.setInternalenergy(eps);
  }

  /**
  * @brief      Compute specific internal energy
  *             Uses abar and temperature
  *
  * @param      particle
  */
  static void
  compute_internal_energy_given_t(body & particle) {
    const double abar = particle.getAbar(),
                 T    = particle.getTemperature();
    double eps = T*KBOL/(AMU*abar*(poly_gamma - 1.));
    particle.setInternalenergy(eps);
  }

 /**
  * @brief      Compute adiabatic invariant: reuse the function from
  *             polytropic EOS
  *
  * @param      particle
  */
  static void
  compute_adiabatic(body & particle){
    eos_t<param::eos_polytropic>::compute_adiabatic(particle);
  }
}; // ...<eos_ideal>


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

  // constants of the piecewise-polytrope fit to the pressure function
  static constexpr double ppt_x0 = 1.25;
  static constexpr double ppt_A1 = 1.6;
  static constexpr double ppt_A2 = 2.0;

public:
  static void
  init(body & particle){
    if(initialize_u) compute_internal_energy(particle);
    particle.setTemperature(0.0);
  }

  static void read_data(){}

  static inline double
  pressure_given_rhoYe(double rho, double Ye) {
    double x = cbrt(rho*Ye/B_wd_nm);
    double x2 = square(x);
    return A_wd*(x*(2.*x2 - 3.)*sqrt(x2 + 1.) + 3.*asinh(x));
  }

  static inline double
  soundspeed_given_rhoYe(double rho, double Ye) {
    double x = cbrt(rho*Ye/B_wd_nm);
    double x2 = square(x);
    double numer = (1. + x2)*(6.*x2 - 3.) + 3. + x2*(2.*x2 - 3.);
    double denom = (1. + x2)*(6.*x2 + 1.) - 1. + x2*(2.*x2 + 1.);
    return sqrt(numer/(3.*denom)) * C_LIGHT_CGS;
  }

  static void
  compute_pressure(body& particle){
    double rho = particle.getDensity();
    double Ye  = particle.getZbar()/particle.getAbar(); // TODO: shouldn't we be using electron fraction here?
    double P = pressure_given_rhoYe(rho, Ye);
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
    double cs = soundspeed_given_rhoYe(rho, Ye);

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
  * @brief      This EOS is temperature-independent, so the function is empty
  *
  * @param      particle
  */
  static void
  compute_temperature(body & particle) {}

  /**
  * @brief      Compute specific internal energy
  *             Uses piecewise-polytrope approximation
  *
  * @param      particle
  */
  static void
  compute_internal_energy(body & particle) {
    const double
        rho = particle.getDensity(),
        abar = particle.getAbar(),
        zbar = particle.getZbar(),
        Ye = zbar/abar;
    const double x   = cbrt(rho*Ye/B_wd_nm),
                 x2  = square(x),
                 x3  = cube(x);
    const double eps = A_wd/rho*(8.*x3*(sqrt(x2 + 1.) - 1.)
                 - (x*(2.*x2 - 3.)*sqrt(x2 + 1.) + 3.*asinh(x)));
    particle.setInternalenergy(eps);
  }

  /**
  * @brief      Compute specific internal energy
  *             Uses piecewise-polytrope approximation
  *
  * @param      particle
  */
  /*static void
  compute_internal_energy(body & particle) {
    const double rho = particle.getDensity(),
                 Ye  = particle.getZbar()/particle.getAbar();
    const double x02 = ppt_x0*ppt_x0;
    const double delta_eps = ppt_A1*x02*x02*0.25 - ppt_A2*ppt_x0*x02/3.0;
    double x = cbrt(rho*Ye/B_wd_nm);
    double x2 = square(x), x4 = square(x2);
    double eps = 0.0;
    if (x < ppt_x0) {
      eps = ppt_A1*x4*0.25;
    }
    else {
      eps = ppt_A2*x*x2/3.0 + delta_eps;
    }
    eps*= A_wd;
    particle.setInternalenergy(eps);
  }*/
}; // ...<eos_wd>

template<>
class eos_t<param::eos_ppt>{
  static double rho_thr;    // density threshold

public:
  /**
  * @brief      Compute adiabatic invariant from density and pressure
  *             In the piecewise-polytropic EOS, adiabatic invariant
  *             corresponds to the first polytropic segment K1:
  *
  *              P(rho) = K1\rho^\Gamma1 + K2\rho^\Gamma2
  *
  * @param      particle
  */
  static void
  compute_adiabatic(body & particle){
    eos_t<param::eos_ppt>::rho_thr = param::ppt_density_thr;
    const double rho = particle.getDensity(),
                 P   = particle.getPressure();
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

  /**
  * @brief      Initialized adiabatic invariant (K1)
  *
  * @param      particle
  */
  static void init(body & particle) {
    compute_adiabatic(particle);
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
    double soundspeed = 0.;
    if (rho < rho_thr) {
      soundspeed = sqrt(K1*poly_gamma*pow(rho,poly_gamma - 2.));
    }
    else {
      double K2 = K1*pow(rho_thr, poly_gamma - poly_gamma2);
      soundspeed = sqrt(K2*poly_gamma2*pow(rho,poly_gamma2 - 2.));
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
    double eps = 0.;
    if (rho < rho_thr) {
      eps = K1*pow(rho, poly_gamma - 1.)/(poly_gamma - 1.);
    }
    else {
      double K2 = K1*pow(rho_thr, poly_gamma - poly_gamma2);
      eps = K2*pow(rho,     poly_gamma2 - 1.)/(poly_gamma2 - 1.)
          - K2*pow(rho_thr, poly_gamma2 - 1.)/(poly_gamma2 - 1.)
          + K1*pow(rho_thr, poly_gamma  - 1.)/(poly_gamma  - 1.);
    }
    particle.setInternalenergy(eps);
  }

};

// declare static member of a templated class
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
  // pressure function constants
  static constexpr double A_wd = 6.00288e22;
  static constexpr double B_wd_nm = 9.81011e5;
public:
  static void init(body & particle){
    particle.setTemperature(initial_temp);
    compute_internal_energy(particle);
  }
  static void read_data(){}

  static void compute_pressure(body& particle){
    const double
      rho = particle.getDensity(),
      Ye  = particle.getZbar()/particle.getAbar();
    double u_gas = get_internal_energy_idealgas(particle);
    double P = (poly_gamma - 1.)*rho*u_gas
             + eos_t<param::eos_wd>::pressure_given_rhoYe(rho, Ye);
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
        abar = particle.getAbar(), zbar = particle.getZbar(),
        Ye = zbar/abar, mu = abar*AMU/(zbar + 1.);
    double u_gas = get_internal_energy_idealgas(particle);
    double cs2 = poly_gamma*(poly_gamma - 1.)*u_gas
               + square(eos_t<param::eos_wd>::soundspeed_given_rhoYe(rho,Ye));
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
        rho  = particle.getDensity(),
        abar = particle.getAbar(),
        zbar = particle.getZbar(),
        mu   = abar*AMU/(zbar + 1.),
        u_gas = get_internal_energy_idealgas(particle);
    particle.setTemperature((poly_gamma - 1.)*u_gas*mu/KBOL);
  }

  static void
  compute_internal_energy(body & particle) {
    // TODO: check
    const double rho = particle.getDensity(),
                zbar = particle.getZbar(), abar = particle.getAbar(),
                 Ye  = zbar/abar,          temp = particle.getTemperature(),
                  mu = abar*AMU/(zbar + 1.);
    const double x  = cbrt(rho*Ye/B_wd_nm),
                 x2 = square(x),
                 x3 = cube(x);
    // calculate degenerate int. energy
    const double u_deg =  A_wd/rho*(8.*x3*(sqrt(x2 + 1.) - 1.)
                 - (x*(2.*x2 - 3.)*sqrt(x2 + 1.) + 3.*asinh(x)));
    // calculate gas int. energy
    const double u_gas = KBOL*temp/mu/(poly_gamma - 1.);
    particle.setInternalenergy(u_deg+u_gas);
  }
private:
  /**
  * @brief      Extracts the ideal gas internal energy from the
  *             int. e value of the particle
  *
  * @param      particle
  */
  static double
  get_internal_energy_idealgas(body & particle) {
    const double
        u   = particle.getInternalenergy(), rho = particle.getDensity(),
        abar = particle.getAbar(), zbar = particle.getZbar(),
        Ye = zbar/abar;
    const double x  = cbrt(rho*Ye/B_wd_nm),
                 x2 = square(x),
                 x3 = cube(x);
    double u_deg = A_wd/rho*(8.*x3*(sqrt(x2 + 1.) - 1.)
                 - (x*(2.*x2 - 3.)*sqrt(x2 + 1.) + 3.*asinh(x)));
    double u_gas = u - u_deg;
    if (u_gas < 0.) u_gas = 0.;
    return u_gas;
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

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

// Root finder (originated from SC reader)
#include "root_finder.h"

// Tabulated EOS utilities and implementations
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

// main eos function type
typedef void (*compute_quantity_t)(body &);

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
  * @brief      Compute "entropy" (actually, adiabatic invariant which is
  *             a function of entropy), from density and pressure
  *
  * @param      particle
  */
  static void
  compute_entropy(body & particle){
    const double rho = particle.getDensity(),
                 P = particle.getPressure();
    double K = adiabatic_given_rhoP(rho, P);
    particle.setEntropy(K);
  }

  /**
  * @brief      Initialize equation of state (nothing for this eos type)
  */
  static void init() {}

  /**
  * @brief      Compute pressure from density using polytrope
  *             P(\rho) = A*\rho^\Gamma
  *
  * @param      particle
  */
  static void compute_pressure(body & particle) {
    const double rho = particle.getDensity(),
                 K   = particle.getEntropy();
    particle.setPressure(K*pow(rho, poly_gamma));
  }

  /**
  * @brief      Compute sound speed for ideal fluid or polytropic eos
  *             cs = sqrt{ A*\Gamma\rho^(\Gamma-1) }
  *
  * @param      particle
  */
  static void compute_soundspeed(body & particle) {
    const double rho = particle.getDensity(),
                 K   = particle.getEntropy();
    double soundspeed = sqrt(K*poly_gamma*pow(rho, poly_gamma - 1.));
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
                 K   = particle.getEntropy();
    double eps = K*pow(rho, poly_gamma - 1.)/(poly_gamma - 1.);
    particle.setInternalenergy(eps);
  }

  compute_quantity_t compute_spct_given_rho_u = nullptr;

}; // ...<eos_polytropic>


template<>
class eos_t<param::eos_ideal>{

public:
  /**
  * @brief      Initialize equation of state (nothing for this eos type)
  */
  static void init() {}

  /**
  * @brief      Computes pressure using the density and internal energy
  *
  * @param      particle
  */
  static void
  compute_pressure(body& particle){
    double pressure =
      (poly_gamma - 1)*particle.getDensity()*particle.getInternalenergy();
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
                 K   = particle.getEntropy();
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
  * @brief      Compute adiabatic invariant (a function of entropy):
  *             reuse the function from polytropic EOS
  *
  * @param      particle
  */
  static void
  compute_entropy(body & particle){
    eos_t<param::eos_polytropic>::compute_entropy(particle);
  }

  compute_quantity_t compute_spct_given_rho_u = nullptr;

}; // ...<eos_ideal>


/**
* @brief      Equation of state for a cold white dwarf.
*             See Chandrasekhar 1939, Ch.11, or
*             Arnett 1996, "Supernovae and Nucleosynthesis", (B.33-34)
*/
template<>
class eos_t<param::eos_wd>{

public:
  // pressure function constants
  // here \lambda_e := h/(m_e c) -- de Broglie wavelength of an electron
  static constexpr double
    A_wd    = 6.00233181e22, // [dynes/cm^2] A_wd = pi/3 m_e c^2/\lambda_e^3
    B_wd_nm = 9.73932099e5;  // [moles/cm^3] B_wd = 8pi / (3 N_A \lambda_e^3)

  // constants of the piecewise-polytrope fit to the pressure function:
  //
  //   P(x)/A_wd = x*(2x^2 - 3)*sqrt(1 + x^2) + 3*asinh(x)
  //
  // which can be approximated by the following:
  //
  //   P(x)/A_wd ~ (x<1.25) ? (1.6*x**5) : (2.0*x**4)
  //
  // This is a piecewise polytrope with parameters:
  //
  //   ppt_density_thr = B_wd/Y_e * ppt_x0**3 = 3.80442e+06*(0.5/Ye) [g/cm3]
  //   poly_gamma  = 5/3
  //   poly_gamma2 = 4/3
  //
  static constexpr double ppt_x0 = 1.25;
  static constexpr double ppt_A1 = 1.6;
  static constexpr double ppt_A2 = 2.0;

  /**
  * @brief      Initialize equation of state (nothing for this eos type)
  */
  static void init() {}

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
    return sqrt(8.*A_wd*Ye*x2/3./B_wd_nm/sqrt(x2 + 1.));
  }

  static void
  compute_pressure(body& particle){
    double rho = particle.getDensity();
    double Ye  = particle.getElectronfraction();
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
    double Ye  = particle.getElectronfraction();
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
  * @brief      Compute entropy
  *             TODO: implement
  *
  * @param      particle
  */
  static void
  compute_entropy(body & particle) {
    /* ... */
  }

  /**
  * @brief      This EOS is temperature-independent, so the function is empty
  *
  * @param      particle
  */
  static void
  compute_temperature(body & particle) {}

  /**
  * @brief      Compute specific internal energy
  *
  * @param      particle
  */
  static void
  compute_internal_energy(body & particle) {
    const double
        rho = particle.getDensity(),
        Ye = particle.getElectronfraction();
    const double x   = cbrt(rho*Ye/B_wd_nm),
                 x2  = square(x),
                 x3  = cube(x);
    const double eps = A_wd/rho*(8.*x3*(sqrt(x2 + 1.) - 1.)
                 - (x*(2.*x2 - 3.)*sqrt(x2 + 1.) + 3.*asinh(x)));
    particle.setInternalenergy(eps);
  }

  compute_quantity_t compute_spct_given_rho_u = nullptr;

}; // ...<eos_wd>

template<>
class eos_t<param::eos_ppt>{
  static double rho_thr;    // density threshold

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
    particle.setEntropy(K1);
  }

  /**
  * @brief      Initialize equation of state (nothing for this eos type)
  */
  static void init() {}

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
                 K1  = particle.getEntropy(),
                 gam = (rho < rho_thr ? poly_gamma : poly_gamma2);
    double soundspeed = 0.;
    if (rho < rho_thr) {
      soundspeed = sqrt(K1*poly_gamma*pow(rho,poly_gamma - 1.));
    }
    else {
      double K2 = K1*pow(rho_thr, poly_gamma - poly_gamma2);
      soundspeed = sqrt(K2*poly_gamma2*pow(rho,poly_gamma2 - 1.));
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
    else {
      double K2 = K1*pow(rho_thr, poly_gamma - poly_gamma2);
      eps = K2*pow(rho,     poly_gamma2 - 1.)/(poly_gamma2 - 1.)
          - K2*pow(rho_thr, poly_gamma2 - 1.)/(poly_gamma2 - 1.)
          + K1*pow(rho_thr, poly_gamma  - 1.)/(poly_gamma  - 1.);
    }
    particle.setInternalenergy(eps);
  }

  compute_quantity_t compute_spct_given_rho_u = nullptr;

};

// declare static member of a templated class
template<>
double eos_t<param::eos_ppt>::rho_thr;

template<>
class eos_t<param::eos_no_eos>{
public:
  static void init(){}
  static void compute_pressure(body& particle){}
  static void compute_soundspeed(body& particle){}
  static void compute_entropy(body& particle){}
  static void compute_temperature(body& particle){}
  static void compute_internal_energy(body& particle){}
};

template<>
class eos_t<param::eos_wd_ideal_gas>{
  // pressure function constants
  static constexpr double A_wd = eos_t<param::eos_wd>::A_wd;
  static constexpr double B_wd_nm = eos_t<param::eos_wd>::B_wd_nm;
public:
  /**
  * @brief      Initialize equation of state (nothing for this eos type)
  */
  static void init() {}

  static void compute_pressure(body& particle){
    const double
      rho = particle.getDensity(),
      Ye  = particle.getElectronfraction();
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
        Ye = particle.getElectronfraction(),
        abar = particle.getAbar(),
        zbar = abar*Ye, mu = abar*AMU/(zbar + 1.);
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
        zbar = abar*particle.getElectronfraction(),
        mu   = abar*AMU/(zbar + 1.),
        u_gas = get_internal_energy_idealgas(particle);
    particle.setTemperature((poly_gamma - 1.)*u_gas*mu/KBOL);
  }

  static void
  compute_internal_energy(body & particle) {
    // TODO: check
    const double rho = particle.getDensity(),
                abar = particle.getAbar(),
                  Ye = particle.getElectronfraction(),
                zbar = abar*Ye,
                temp = particle.getTemperature(),
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

  compute_quantity_t compute_spct_given_rho_u = nullptr;

private:
  /**
  * @brief      Extracts the ideal gas internal energy from the
  *             int. e value of the particle
  *
  * @param      particle
  */
  static double
  get_internal_energy_idealgas(body & particle) {
    const double rho = particle.getDensity(),
                   u = particle.getInternalenergy(),
                abar = particle.getAbar(),
                  Ye = particle.getElectronfraction(),
                zbar = abar*Ye;
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

#ifdef eos_type
#  define compute_pressure     eos_t<eos_type>::compute_pressure
#  define compute_soundspeed   eos_t<eos_type>::compute_soundspeed
#  define compute_entropy      eos_t<eos_type>::compute_entropy
#  define compute_temperature  eos_t<eos_type>::compute_temperature
#  define compute_internal_energy eos_t<eos_type>::compute_internal_energy
#  define compute_spct_given_rho_u eos_t<eos_type>::compute_spct_given_rho_u
#else
compute_quantity_t compute_pressure = nullptr;
compute_quantity_t compute_soundspeed = nullptr;
compute_quantity_t compute_entropy = nullptr;
compute_quantity_t compute_temperature = nullptr;
compute_quantity_t compute_internal_energy = nullptr;
compute_quantity_t compute_spct_given_rho_u = nullptr;
#endif

/**
 * @brief  Installs the 'compute_pressure' and 'compute_soundspeed'
 *         function pointers, depending on the value of eos_type
 */
void
select() {
  using namespace param;

#ifndef eos_type
  log_one(info) << "Selecting equation of state: "
                << eos_type_decode[eos_type] << std::endl;
  switch(eos_type){
    case(eos_ideal):
      compute_pressure = eos_t<eos_ideal>::compute_pressure;
      compute_soundspeed = eos_t<eos_ideal>::compute_soundspeed;
      compute_entropy = eos_t<eos_ideal>::compute_entropy;
      compute_temperature = eos_t<eos_ideal>::compute_temperature;
      compute_internal_energy = eos_t<eos_ideal>::compute_internal_energy;
      eos_t<eos_ideal>::init();
      break;
    case(eos_polytropic):
      compute_pressure = eos_t<eos_polytropic>::compute_pressure;
      compute_soundspeed = eos_t<eos_polytropic>::compute_soundspeed;
      compute_entropy = eos_t<eos_polytropic>::compute_entropy;
      compute_temperature = eos_t<eos_polytropic>::compute_temperature;
      compute_internal_energy = eos_t<eos_polytropic>::compute_internal_energy;
      eos_t<eos_polytropic>::init();
      break;
    case(eos_wd):
      compute_pressure = eos_t<eos_wd>::compute_pressure;
      compute_soundspeed = eos_t<eos_wd>::compute_soundspeed;
      compute_entropy = eos_t<eos_wd>::compute_entropy;
      compute_temperature = eos_t<eos_wd>::compute_temperature;
      compute_internal_energy = eos_t<eos_wd>::compute_internal_energy;
      eos_t<eos_wd>::init();
      break;
    case(eos_ppt):
      compute_pressure = eos_t<eos_ppt>::compute_pressure;
      compute_soundspeed = eos_t<eos_ppt>::compute_soundspeed;
      compute_entropy = eos_t<eos_ppt>::compute_entropy;
      compute_temperature = eos_t<eos_ppt>::compute_temperature;
      compute_internal_energy = eos_t<eos_ppt>::compute_internal_energy;
      eos_t<eos_ppt>::init();
      break;
    case(eos_no_eos):
      compute_pressure = eos_t<eos_no_eos>::compute_pressure;
      compute_soundspeed = eos_t<eos_no_eos>::compute_soundspeed;
      compute_entropy = eos_t<eos_no_eos>::compute_entropy;
      compute_temperature = eos_t<eos_no_eos>::compute_temperature;
      compute_internal_energy = eos_t<eos_no_eos>::compute_internal_energy;
      break;
    case(eos_stellar_collapse):
      compute_pressure = eos_t<eos_stellar_collapse>::compute_pressure;
      compute_soundspeed = eos_t<eos_stellar_collapse>::compute_soundspeed;
      compute_entropy = eos_t<eos_stellar_collapse>::compute_entropy;
      compute_temperature = eos_t<eos_stellar_collapse>::compute_temperature;
      compute_internal_energy = eos_t<eos_stellar_collapse>::compute_internal_energy;
      eos_t<eos_stellar_collapse>::init();
      break;
    case(eos_wd_ideal_gas):
      compute_pressure = eos_t<eos_wd_ideal_gas>::compute_pressure;
      compute_soundspeed = eos_t<eos_wd_ideal_gas>::compute_soundspeed;
      compute_temperature = eos_t<eos_wd_ideal_gas>::compute_temperature;
      compute_internal_energy = eos_t<eos_wd_ideal_gas>::compute_internal_energy;
      eos_t<eos_wd_ideal_gas>::init();
      break;
    case(eos_helmholtz):
      compute_pressure = eos_t<eos_helmholtz>::compute_pressure;
      compute_soundspeed = eos_t<eos_helmholtz>::compute_soundspeed;
      compute_entropy = eos_t<eos_helmholtz>::compute_entropy;
      compute_temperature = eos_t<eos_helmholtz>::compute_temperature;
      compute_internal_energy = eos_t<eos_helmholtz>::compute_internal_energy;
      eos_t<eos_helmholtz>::init();

      //// Check Helmholtz table and root finder
      //// eos_t<eos_helmholtz>::table_check();
      //eos_t<eos_helmholtz>::root_finder_check();
      //MPI_Finalize();
      //exit(0);

      break;
    default:
      std::cerr<<"Undefined eos type"<<std::endl;
      MPI_Finalize();
      exit(0);
  }
#endif // eos_type
} // select

/**
 * @brief      set uniform average atomic weight (abar) and electron
 *             fraction Ye := zbar/abar, using parameters initial_abar and
 *             initial_zbar
 *             TODO: read from species file
 *
 * @param      particle
 *
 * @uses       initial_abar     global parameter
 * @uses       initial_zbar     global parameter
 */
void
initialize_abarzbar(body & particle)
{
  using namespace param;
  particle.setAbar(initial_abar);
  particle.setElectronfraction(initial_zbar/initial_abar);
} // initialize_abarzbar

} // namespace eos

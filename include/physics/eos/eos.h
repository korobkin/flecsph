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
#include "units.h"
#include "eos_utils.h"
#include "eos_ppt.h"
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
  * @brief      Initialize equation of state (nothing for this eos type)
  */
  static void init() {}

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

  // TODO
  static get_quantity_t get_dpdrho_at_temp;
  static compute_quantity_t compute_spct_given_rho_u;

}; // ...<eos_polytropic>

#if eos_type == eos_polytropic
  get_quantity_t eos_t<param::eos_polytropic>::get_dpdrho_at_temp = nullptr;
  compute_quantity_t eos_t<param::eos_polytropic>::compute_spct_given_rho_u = nullptr;
#endif


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
    double T = phys::amu/phys::kB*abar*(poly_gamma - 1.)*eps;
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
    double eps = T*phys::kB/(phys::amu*abar*(poly_gamma - 1.));
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

  // TODO
  static get_quantity_t get_dpdrho_at_temp;
  static compute_quantity_t compute_spct_given_rho_u;

}; // ...<eos_ideal>

#if eos_type == eos_ideal
  get_quantity_t eos_t<param::eos_ideal>::get_dpdrho_at_temp = nullptr;
  compute_quantity_t eos_t<param::eos_ideal>::compute_spct_given_rho_u = nullptr;
#endif


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
  static constexpr double ppt_A1 = 1.5999775;
  static constexpr double ppt_A2 = 2.0;

  /**
  * @brief      Initialize equation of state
  */
  static void init() {}

  static inline double
  pressure_given_rhoYe(double rho, double Ye) {
    double x = cbrt(rho*Ye/B_wd_nm);
    double x2 = square(x);
    return A_wd*((x>0.01) ? (x*(2.*x2 - 3.)*sqrt(x2 + 1.) + 3.*asinh(x))
               : (ppt_A1*x*x2*x2)); // accurate asymptotic at low x
  }

  static inline double
  dPdrho_given_rhoYe(double rho, double Ye) {
    double x = cbrt(rho*Ye/B_wd_nm);
    double x2 = square(x);
    return 8.*A_wd*Ye*x2/3./B_wd_nm/sqrt(x2 + 1.);
  }

  static inline double
  eint_given_rhoYe(double rho, double Ye) {
    double x3  = rho*Ye/B_wd_nm,
           x   = cbrt(x3),
           x2  = x*x;
    return A_wd/rho*((x>0.01) ? (8.*x3*(sqrt(x2 + 1.) - 1.)
               - (x*(2.*x2 - 3.)*sqrt(x2 + 1.) + 3.*asinh(x)))
               : (1.5*ppt_A1*x3*x2));
  }

  static inline double
  soundspeed_given_rhoYe(double rho, double Ye) {
    return sqrt(dPdrho_given_rhoYe(rho, Ye));
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
    double rho = particle.getDensity();
    double Ye  = particle.getElectronfraction();
    double eps = eint_given_rhoYe(rho, Ye);
    particle.setInternalenergy(eps);
  }

  // TODO
  static get_quantity_t get_dpdrho_at_temp;
  static compute_quantity_t compute_spct_given_rho_u;

}; // ...<eos_wd>

#if eos_type == eos_wd
  get_quantity_t eos_t<param::eos_wd>::get_dpdrho_at_temp = nullptr;
  compute_quantity_t eos_t<param::eos_wd>::compute_spct_given_rho_u = nullptr;
#endif


//*****************************************************************************************
template<>
class eos_t<param::eos_no_eos>{
public:
  static void init() {}
  static void compute_pressure(body& particle){}
  static void compute_soundspeed(body& particle){}
  static void compute_entropy(body& particle){}
  static void compute_temperature(body& particle){}
  static void compute_internal_energy(body& particle){}
  static get_quantity_t get_dpdrho_at_temp;
  static compute_quantity_t compute_spct_given_rho_u;
};
#if eos_type == eos_no_eos
  get_quantity_t eos_t<param::eos_no_eos>::get_dpdrho_at_temp = nullptr;
  compute_quantity_t eos_t<param::eos_no_eos>::compute_spct_given_rho_u = nullptr;
#endif

template<>
class eos_t<param::eos_wd_thermal>{
  // pressure function constants
  static constexpr double A_wd = eos_t<param::eos_wd>::A_wd;
  static constexpr double B_wd_nm = eos_t<param::eos_wd>::B_wd_nm;
public:
  /**
  * @brief      Initialize equation of state
  */
  static void init() {}

  static inline double
  eint_given_rho_temp(const double rho, const double temp,
      const double abar, const double zbar) {
    double temp2 = temp*temp;
    double u_ph = phys::arad*temp2*temp2/rho;
    double u_ions = 1.5*phys::Rgas*(zbar + 1.)/abar * temp;
    double u_deg = eos_t<param::eos_wd>::eint_given_rhoYe(rho, zbar/abar);
    return u_ph + u_ions + u_deg;
  }

  static inline double
  pressure_given_rho_temp(const double rho, const double temp,
      const double abar, const double zbar) {
    double temp2 = temp*temp;
    double P_ph = (1./3.)*phys::arad*temp2*temp2;
    double P_ions = phys::Rgas*(zbar + 1.)/abar * rho*temp;
    double P_deg = eos_t<param::eos_wd>::pressure_given_rhoYe(rho, zbar/abar);
    return P_ph + P_ions + P_deg;
  }

  static double
  temp_given_rho_eint(const double rho, const double eint,
      const double abar, const double zbar) {
    const double temperature_floor = 1000.;
    double temp = temperature_floor;

    // subtract degenerate energy (only depends on density)
    double u = eint - eos_t<param::eos_wd>::eint_given_rhoYe(rho, zbar/abar);
    if (u > 0.) {
      // initial guess
      temp = sqrt(sqrt(rho*u/phys::arad));

      // a few newton-raphsons
      for (int i = 0; i < 5; ++i) {
        double temp2 = temp*temp;
        double du1dT = 4.*phys::arad*temp*temp2/rho;
        double du2dT = 1.5*phys::Rgas*(zbar + 1.)/abar;
        double eint = (0.25*du1dT + du2dT)*temp - u;
        double delta_temp = eint/(du1dT + du2dT);
        if (std::abs(delta_temp/temp) < 1e-12)
          break;
        temp -= delta_temp;
      }
    } // else (if residual u < 0), return temperature floor
    return temp;
  }

  static double
  pressure_given_rho_eint(const double rho, const double eint,
      const double abar, const double zbar) {
    double temp = temp_given_rho_eint(rho, eint, abar, zbar);
    return    pressure_given_rho_temp(rho, temp, abar, zbar);
  }

  static double
  soundspeed_given_rho_temp(const double rho, const double T,
      const double abar, const double zbar) {

    double T3 = T*T*T;
    double P = pressure_given_rho_temp(rho, T, abar, zbar);
    double Ye = zbar/abar;
    double dPdd = eos_t<param::eos_wd>::dPdrho_given_rhoYe(rho, Ye)
                + phys::Rgas*(zbar + 1.)/abar * T;
    double dPdT = 4./3.*phys::arad*T3
                + phys::Rgas*(zbar + 1.)/abar * rho;
    double dPdT2= dPdT*dPdT;
    double dudT = 4.*phys::arad*T3/rho
                + 1.5*phys::Rgas*(zbar + 1.)/abar;
    double eint = eint_given_rho_temp(rho, T, abar, zbar);
    double denom = 1. + (eint + P/rho)/(phys::clight*phys::clight);
    double numer = dPdd + T/(rho*rho)*dPdT2/dudT;

    return sqrt(numer/denom);
  }

  static double
  soundspeed_given_rho_eint(const double rho, const double eint,
      const double abar, const double zbar) {

    double T = temp_given_rho_eint(rho, eint, abar, zbar);
    return  soundspeed_given_rho_temp(rho, T, abar, zbar);
  }

  /**
  * @brief      Compute pressure given density and internal energy
  *
  * @param      particle
  */
  static void
  compute_pressure(body & particle) {
    const double
      rho = particle.getDensity(),
      eint = particle.getInternalenergy(),
      Ye  = particle.getElectronfraction(),
      abar = particle.getAbar(),
      zbar = abar*Ye;
    double P = pressure_given_rho_eint(rho, eint, abar, zbar);
    // DEBUG
    // if (P != P) {
    //   double T = temp_given_rho_eint(rho, eint, abar, zbar);
    //   std::cout << "ERROR: pressure is NaN" << std::endl;
    //   std::cout << "Failed particle id: " << particle.id() << std::endl;
    //   std::cerr << "particle position: " << particle.coordinates() << std::endl;
    //   std::cerr << "particle velocity: " << particle.getVelocity() << std::endl;
    //   std::cerr << "particle acceleration: " << particle.getAcceleration() << std::endl;
    //   std::cerr << "particle density: " << particle.getDensity() << std::endl;
    //   std::cerr << "particle internal energy: " << particle.getInternalenergy() << std::endl;
    //   std::cerr << "particle expected temperature: " << particle.getTemperature() << std::endl;
    //   std::cerr << "particle computed temperature: " << T << std::endl;
    //   std::cerr << "particle ye: " << particle.getElectronfraction() << std::endl;
    //   std::cerr << "particle abar: " << particle.getAbar() << std::endl;
    //   std::cerr << "smoothing length:  " << particle.radius() << std::endl;
    //   MPI_Abort(MPI_COMM_WORLD, -1);
    // }
    particle.setPressure(P);
  }

  /**
  * @brief      Compute soundspeed from density and internal energy
  *
  * @param      particle
  */
  static void
  compute_soundspeed(body & particle) {
    const double
        rho = particle.getDensity(),
        eint = particle.getInternalenergy(),
        Ye = particle.getElectronfraction(),
        abar = particle.getAbar(),
        zbar = abar*Ye;

    double cs = soundspeed_given_rho_eint(rho, eint, abar, zbar);
    particle.setSoundspeed(cs);
  }

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
  * @brief      Compute temperature from density and internal energy
  *
  * @param      particle
  */
  static void
  compute_temperature(body & particle) {
    const double
        rho = particle.getDensity(),
        eint = particle.getInternalenergy(),
        Ye = particle.getElectronfraction(),
        abar = particle.getAbar(),
        zbar = abar*Ye;
    double T = temp_given_rho_eint(rho, eint, abar, zbar);
    particle.setTemperature(T);
  }

  /**
  * @brief      Compute internal energy given density and temperature
  *
  * @param      particle
  */
  static void
  compute_internal_energy(body & particle) {
    const double rho = particle.getDensity(),
                abar = particle.getAbar(),
                  Ye = particle.getElectronfraction(),
                zbar = abar*Ye,
                temp = particle.getTemperature();
    double eint = eint_given_rho_temp(rho, temp, abar, zbar);
    particle.setInternalenergy(eint);
  }

  /**
  * @brief      Returns (dP/drho)_T: partial derivatie of the pressure
  *             with respect to density at constant temperature
  *
  * @param      particle
  */
  static inline double
  get_dpdrho_at_temp(const body & particle) {
    const double
      rho = particle.getDensity(),
      temp= particle.getTemperature(),
      Ye  = particle.getElectronfraction(),
      abar = particle.getAbar(),
      zbar = abar*Ye;
    double dP_ph = 0.;
    double dP_ions = phys::Rgas*(zbar + 1.)/abar * temp;
    double dP_deg = eos_t<param::eos_wd>::dPdrho_given_rhoYe(rho, Ye);
    return dP_ph + dP_ions + dP_deg;
  }

  static compute_quantity_t compute_spct_given_rho_u;

}; // ...<eos_wd_thermal>

#if eos_type == eos_wd_thermal
  compute_quantity_t eos_t<param::eos_wd_thermal>::compute_spct_given_rho_u = nullptr;
#endif

#ifdef eos_type
#  define compute_pressure     eos_t<eos_type>::compute_pressure
#  define compute_soundspeed   eos_t<eos_type>::compute_soundspeed
#  define compute_entropy      eos_t<eos_type>::compute_entropy
#  define compute_temperature  eos_t<eos_type>::compute_temperature
#  define compute_internal_energy eos_t<eos_type>::compute_internal_energy
#  define compute_spct_given_rho_u eos_t<eos_type>::compute_spct_given_rho_u
#  define get_dpdrho_at_temp   eos_t<eos_type>::get_dpdrho_at_temp
#else
compute_quantity_t compute_pressure = nullptr;
compute_quantity_t compute_soundspeed = nullptr;
compute_quantity_t compute_entropy = nullptr;
compute_quantity_t compute_temperature = nullptr;
compute_quantity_t compute_internal_energy = nullptr;
compute_quantity_t compute_spct_given_rho_u = nullptr;
get_quantity_t get_dpdrho_at_temp = nullptr;
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
      eos_t<eos_no_eos>::init();
      break;
    case(eos_stellar_collapse):
      compute_pressure = eos_t<eos_stellar_collapse>::compute_pressure;
      compute_soundspeed = eos_t<eos_stellar_collapse>::compute_soundspeed;
      compute_entropy = eos_t<eos_stellar_collapse>::compute_entropy;
      compute_temperature = eos_t<eos_stellar_collapse>::compute_temperature;
      compute_internal_energy = eos_t<eos_stellar_collapse>::compute_internal_energy;
      compute_spct_given_rho_u = eos_t<eos_stellar_collapse>::compute_spct_given_rho_u;
      eos_t<eos_stellar_collapse>::init();
      break;
    case(eos_wd_thermal):
      compute_pressure = eos_t<eos_wd_thermal>::compute_pressure;
      compute_soundspeed = eos_t<eos_wd_thermal>::compute_soundspeed;
      compute_entropy = eos_t<eos_wd_thermal>::compute_entropy;
      compute_temperature = eos_t<eos_wd_thermal>::compute_temperature;
      compute_internal_energy = eos_t<eos_wd_thermal>::compute_internal_energy;
      get_dpdrho_at_temp = eos_t<eos_wd_thermal>::get_dpdrho_at_temp;
      eos_t<eos_wd_thermal>::init();
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

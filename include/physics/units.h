/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2017 Triad National Security, LLC
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
 * @file units.h
 * @author Oleg Korobkin
 * @date June 2019
 * @brief Handle simulation units and physical constants
 *
 * In a multiphysics simulation, it is important to keep a record of units and
 * make sure they are used uniformly across different modules. In FleCSPH,
 * default units are CGS, but other systems can be set using the parameter
 *
 *   simulation_units = "[cgs|si|geom|nuc]"
 *
 * The options (implemented so far) are as follows:
 *  - cgs:  the CGS system;
 *  - si:   SI system (2019 standard)
 *  - geom: geometric units with G=c=1, [M]=Msun;
 *  - nuc:  nuclear units: [T]=fs, [Energy]=MeV, [L]=fm
 */


#ifndef UNITS_H
#define UNITS_H

#define UNITS_H_STR_HELPER(x) #x
#define UNITS_H_STR(x) UNITS_H_STR_HELPER(x)

#include <iostream>
#include <cmath>
#include "params.h"

//////////////////////////////////////////////////////////////////////
namespace phys {

// basic units: what numerical units correspond to in the simulation units;
// for example, if simulation_units = "cgs" and units::length = 2.0, it means
// that numerical length of 1 corresponds to 2cm.
constexpr double 
#if (units == CGS_UNITS)
    M = 1.0, // mass
    L = 1.0, // length
    T = 1.0, // time
    K = 1.0; // temperature
#elif (units == SI_UNITS)
    M = 1000.0, // kg = 1000 g
    L = 100.0,  // m = 100 cm
    T = 1.0,    // time
    K = 1.0;    // temperature
#elif (units == GEOM_UNITS)
    M = 1.989e33,       // Solar mass [g]
    L = 1.476961476e+5, // grav. radius of a Solar-mass BH: GM/c^2 [cm]
    T = 4.92661318e-06, // time to cross grav. radius [s]
    K = 1.0;            // temperature (shall we switch to MeVs?)
#endif
// macro debugging message:
// #pragma message "from units.h: units = " UNITS_H_STR(units)

// derived quantities
constexpr double
    T2 = T*T,
    L2 = L*L,
    L3 = L*L*L,
    K4 = K*K*K*K;

constexpr double
    F = M*L/T2,       // unit of force
    E = M*L2/T2,      // unit of energy
    Q = sqrt(M*L3)/T; // charge

////////
// The following structure is used in order to convert foreign
// tables that are written in other units, to the current
// code units. They provide converstion factors for various physical
// quantities to the current system of units.
// Usage example:
// 
//  mass_ = mass_in_SI * phys::cfactor_to_<param::si_units>::mass
//
template<param::units_keyword_enum U>
struct cfactor_to_ { };

template<>
struct cfactor_to_<param::cgs_units> {
  static constexpr double
      mass = 1./M,
      length = 1./L,
      time = 1./T,
      force = 1./F,
      energy = 1./E,
      velocity = 1./L*T,
      density = 1./M*L3,
      volume = 1./L3,
      temperature = 1./K,
      charge = 1./Q;
};

template<>
struct cfactor_to_<param::si_units> {
  static constexpr double
      mass = 1000./M,
      length = 100./L,
      time = 1./T,
      force = 1e+5/F,
      energy = 1e+7/E,
      velocity = 100./L*T,
      density = 1e-3/M*L3,
      volume = 1e6/L3,
      temperature = 1./K,
      charge = std::sqrt(1e9)/Q;
};

template<>
struct cfactor_to_<param::geom_units> {
  static constexpr double
      mass     = 1.989e+33/M,
      length   = 1.476961476e+5/L,
      time     = 4.92661318e-06/T,
      force    = 1.21033898e+49/F,
      energy   = 1.78762405e+54/E,
      velocity = 2.997924583e+10/L*T,
      density  = 6.173440692e+17/M*L3,
      volume   = 3.221866216e+15/L3,
      temperature = 1./K,
      charge   = 5.138338117e+29/Q;
};

// physics constants: converted from CGS
constexpr double
   clight    = 2.99792458e10  /L*T,     // speed of light
   qe        = 4.80320680e-10 /Q,       // elementary charge
   me        = 9.1093826e-28  /M,       // electron mass
   mp        = 1.67262171e-24 /M,       // proton mass
   mn        = 1.67492728e-24 /M,       // neutron mass
   amu       = 1.66053878e-24 /M,       // Atomic Mass Unit [g/baryon]
   hplanck   = 6.6260693e-27  /E/T,     // Planck constant
   hbar      = 1.0545717e-27  /E/T,     // barred Planck (h over 2pi) 
   kB        = 1.3806505e-16  /E*K,     // Boltzmann constant
   GN        = 6.67384e-8     /L/E*(M*M),   // Newton's gravity constant
   sB        = 5.670400e-5    /E*(L2*T*K4), // Stephan-Boltzmann constant
   arad      = 7.5657e-15     /E*(L3*K4),   // radiation density constant
   thomson   = 6.65245873e-25 /L2,          // Thomson scattering cross-section
   jy        = 1.e-23         /E*(L2*T),    // Jansky (luminous flux)
   pc        = 3.085678e18    /L,       // parsec
   au        = 1.49597870691e13/L,      // a.u. = astronomical unit
   Msun      = 1.989e33       /M,       // solar mass
   Rsun      = 6.96e10        /L,       // solar radius
   Lsun      = 3.827e33       /E*T,     // solar power
   eV        = 1.60217653e-12 /E,       // electronvolt
   MeV       = 1.60217653e-6  /E,       // mega-electronvolt 
   GeV       = 1.60217653e-3  /E,       // giga-electronvolt 
   NAvo      = 6.0221417930e23,         // Avogadro's number [mol^-1]
   Rgas      = kB*NAvo;                 // Ideal gas constant [erg K^-1 mol^-1]

void
output_constants() {
  std::cout << 
    "speed of light                    (clight)  = " << clight  << "\n"
    "elementary charge                 (qe)      = " << qe      << "\n"
    "electron mass                     (me)      = " << me      << "\n"
    "proton mass                       (mp)      = " << mp      << "\n"
    "neutron mass                      (mn)      = " << mn      << "\n"
    "Planck constant                   (hplanck) = " << hplanck << "\n"
    "barred Planck (h over 2pi)        (hbar)    = " << hbar    << "\n"
    "Boltzmann constant                (kB)      = " << kB      << "\n"
    "Newton's gravity constant         (GN)      = " << GN      << "\n"
    "Stephan-Boltzmann constant        (sB)      = " << sB      << "\n"
    "radiation density constant        (arad)    = " << arad    << "\n"
    "Thomson scattering cross-section  (thomson) = " << thomson << "\n"
    "Jansky (luminous flux)            (jy)      = " << jy      << "\n"
    "parsec                            (pc)      = " << pc      << "\n"
    "a.u. = astronomical unit          (au)      = " << au      << "\n"
    "solar mass                        (Msun)    = " << Msun    << "\n"
    "solar radius                      (Rsun)    = " << Rsun    << "\n"
    "solar power                       (Lsun)    = " << Lsun    << "\n"
    "electronvolt                      (eV)      = " << eV      << "\n"
    "mega-electronvolt                 (MeV)     = " << MeV     << "\n"
    "giga-electronvolt                 (GeV)     = " << GeV     << "\n"
  ;
}

template<param::units_keyword_enum U> 
void
output_conversion_factors() {
//  mass_ = mass_in_SI * phys::cfactor_to_<param::si_units>::mass
  std::cout << 
    "conversion factor from " << U << " to current system of units: \n" 
    "               mass *= " << cfactor_to_<U>::mass         << "\n"
    "             length *= " << cfactor_to_<U>::length       << "\n"
    "               time *= " << cfactor_to_<U>::time         << "\n"
    "              force *= " << cfactor_to_<U>::force        << "\n"
    "             energy *= " << cfactor_to_<U>::energy       << "\n"
    "           velocity *= " << cfactor_to_<U>::velocity     << "\n"
    "            density *= " << cfactor_to_<U>::density      << "\n"
    "             volume *= " << cfactor_to_<U>::volume       << "\n"
    "        temperature *= " << cfactor_to_<U>::temperature  << "\n"
    "             charge *= " << cfactor_to_<U>::charge       << "\n"
  ;
}

#if (units == CGS_UNITS)
#  define LENGTH_UNIT_STR   "cm"
#  define MASS_UNIT_STR     "g"
#  define TIME_UNIT_STR     "s"
#  define DENSITY_UNIT_STR  "g/cm3"
#  define PRESSURE_UNIT_STR "dynes/cm2"
#elif (units == SI_UNITS)
#  define LENGTH_UNIT_STR   "m"
#  define MASS_UNIT_STR     "kg"
#  define TIME_UNIT_STR     "s"
#  define DENSITY_UNIT_STR  "kg/m3"
#  define PRESSURE_UNIT_STR "Pa"
#elif (units == GEOM_UNITS)
#  define LENGTH_UNIT_STR   "M"
#  define MASS_UNIT_STR     "M"
#  define TIME_UNIT_STR     "M"
#  define DENSITY_UNIT_STR  "M^-2"
#  define PRESSURE_UNIT_STR "M^-2"
#endif

} // namespace phys
#undef UNITS_H_STR_HELPER(x) #x
#undef UNITS_H_STR(x) UNITS_H_STR_HELPER(x)
#endif // UNITS_H

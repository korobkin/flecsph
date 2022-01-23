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

#include <cmath>
#include "params.h"

//////////////////////////////////////////////////////////////////////
namespace phys {

// basic units: what numerical units correspond to in the simulation units;
// for example, if simulation_units = "cgs" and units::length = 2.0, it means
// that numerical length of 1 corresponds to 2cm.
const double 
#if units == cgs_units
    M = 1.0, // mass
    L = 1.0, // length
    T = 1.0, // time
    K = 1.0; // temperature
#elif units == si_units
    M = 1000.0, // mass
    L = 100.0, // length
    T = 1.0, // time
    K = 1.0; // temperature
#endif

// derived quantities
const double
    T2 = T*T,
    L2 = L*L,
    L3 = L*L*L,
    K4 = K*K*K*K;

const double
    F = M*L/T2,       // unit of force
    E = M*L2/T2,      // unit of energy
    Q = sqrt(M*L3)/T; // charge

// physics constants: default CGS
const double
   clight    = 2.99792458e10  *L/T,     // speed of light
   qe        = 4.80320680e-10 *Q,       // elementary charge
   me        = 9.1093826e-28  *M,       // electron mass
   mp        = 1.67262171e-24 *M,       // proton mass
   mn        = 1.67492728e-24 *M,       // neutron mass
   hplanck   = 6.6260693e-27  *E*T,     // Planck constant
   hbar      = 1.0545717e-27  *E*T,     // barred Planck (h over 2pi) 
   kB        = 1.3806505e-16  *E/K,     // Boltzmann constant
   GN        = 6.67384e-8     *L*E/(M*M),   // Newton's gravity constant
   sB        = 5.670400e-5    *E/(L2*T*K4), // Stephan-Boltzmann constant
   arad      = 7.5657e-15     *E/(L3*K4),   // radiation density constant
   thomson   = 6.65245873e-25 *L2,          // Thomson scattering cross-section
   jy        = 1.e-23         *E/(L2*T),    // Jansky (luminous flux)
   pc        = 3.085678e18    *L,       // parsec
   au        = 1.49597870691e13*L,      // a.u. = astronomical unit
   Msun      = 1.989e33       *M,       // solar mass
   Rsun      = 6.96e10        *L,       // solar radius
   Lsun      = 3.827e33       *E/T,     // solar power
   eV        = 1.60217653e-12 *E,       // electronvolt
   MeV       = 1.60217653e-6  *E,       // mega-electronvolt 
   GeV       = 1.60217653e-3  *E;       // giga-electronvolt 

} // namespace phys
#endif // UNITS_H

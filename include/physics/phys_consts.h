/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2018 Triad National Security, LLC
 * All rights reserved.
 * --------------------------------------------------------------------------~*/

/******************************************************************************
 *                                                                            *
 * PHYS_CONSTS.H                                                              *
 *                                                                            *
 * GLOBAL CONSTANTS AND CONVERSION FACTORS                                    *
 *                                                                            *
 ******************************************************************************/

#pragma once 

#include <math.h>
#include <stdlib.h>

// Fundamental constants in CGS
constexpr double 
  M_SUN_CGS   = 1.98847e33,        // Solar mass [g]
  R_SUN_CGS   = 6.957e10,          // Solar radius [cm]
  C_LIGHT_CGS = 2.99792458e10,     // Speed of light [cm/s]
  EE          = 4.80320680e-10,    // Electron charge [CGS]
  ME          = 9.1093826e-28,     // Electron mass [g]
  MP          = 1.67262171e-24,    // Proton mass [g]
  MN          = 1.67492728e-24,    // Neutron mass [g]
  AMU         = 1.66053878283e-24, // Atomic Mass Unit [g/baryon]
  HPL         = 6.62607015e-27,    // Planck constant [erg*s]
  HBAR        = HPL/(2*M_PI),      // Reduced Planck constant [erg*s]
  KBOL        = 1.3806505e-16,     // Boltzmann constant [erg/K]
  GNEWT       = 6.67259e-8,        // Gravitational constant [cm^3 g^-1 s^-2]
  SIG         = 5.670400e-5,       // Stefan-Boltzmann constant [erg cm^-2 s^-1 K^-4]
  AR          = 4*SIG/C_LIGHT_CGS, // Radiation constant [erg cm^-3 K^-4]
  THOMSON     = 0.665245873e-24,   // Thomson cross section [cm^2]
  COULOMB_LOG = 20.,               // Coulomb logarithm [.]
  ALPHAFS     = 0.007299270073,    // Fine structure constant ~ 1./137. [.]
  GA          = -1.272323,         // Axial-vector coupling [.]
  GA2         = GA*GA,             // Axial-vector coupling squared [.]
  S2THW       = 0.222321,          // sin^2(Theta_W), Theta_W = Weinberg angle 
  S4THW       = S2THW*S2THW,       // sin^4(Theta_W), Theta_W = Weinberg angle
  NUSIGMA0    = 1.7611737037e-44,  // Fundamental neutrino cross section [cm^2]
  AVO         = 6.0221417930e23;   // Avogadro's number [mol^-1]

// nonstandard unit conversion factors
constexpr double 
  EV   = 1.60217653e-12,   // Electron-volt [erg]
  MEV  = 1.0e6 * EV,       // Mega-Electron-Volt [erg]
  GEV  = 1.0e9 * EV,       // Giga-Electron-Volt [erg]
  JY   = 1.e-23,           // Jansky [erg cm^-2 s^-1 Hz^-1]
  PC   = 3.085678e18,      // Parsec [cm]
  AU   = 1.49597870691e13, // Astronomical unit [cm]
  HOUR = 3600.,            // hour [s]
  DAY  = 86400.,           // day [s]
  YEAR = 3.15576e+7;       // Julian year = 365.25 d [s] 

// conversion factors between Geometrized and CGS units 
constexpr double
  C_LIGHT_SQ = C_LIGHT_CGS*C_LIGHT_CGS,         // c^2
  C_LIGHT_QU = C_LIGHT_SQ*C_LIGHT_SQ,           // c^4
  M_GEOM_TO_CGS = M_SUN_CGS,                    // unit of mass
  L_GEOM_TO_CGS = GNEWT*M_SUN_CGS/C_LIGHT_SQ,   // unit of length
  T_GEOM_TO_CGS = L_GEOM_TO_CGS/C_LIGHT_CGS,    // unit of time
  RHO_GEOM_TO_CGS = M_GEOM_TO_CGS               // unit of density
                  /(L_GEOM_TO_CGS*L_GEOM_TO_CGS*L_GEOM_TO_CGS),
  VEL_GEOM_TO_CGS = C_LIGHT_CGS,                // unit of velocity
  P_GEOM_TO_CGS = RHO_GEOM_TO_CGS*C_LIGHT_SQ,   // unit of pressure
  EN_GEOM_TO_CGS = M_GEOM_TO_CGS*C_LIGHT_SQ,    // unit of energy
  EPS_GEOM_TO_CGS = C_LIGHT_SQ,                 // unit of specific energy
  ACC_GEOM_TO_CGS = C_LIGHT_SQ/L_GEOM_TO_CGS;   // unit of acceleration


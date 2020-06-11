/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2018 Triad National Security, LLC
 * All rights reserved.
 * --------------------------------------------------------------------------~*/

/******************************************************************************
 *                                                                            *
 * EOS_CONSTS.H                                                               *
 *                                                                            *
 * GLOBAL CONSTANTS                                                           *
 *                                                                            *
 ******************************************************************************/

#pragma once 

#include <math.h>
#include <stdlib.h>

// Fundamental constants in CGS
constexpr double M_SUN_CGS   = 1.98847e33;            // Solar mass in CGS
constexpr double C_LIGHT_CGS = 2.99792458e10;         // Speed of light in CGS
constexpr double EE          = 4.80320680e-10;        // Electron charge
constexpr double ME          = 9.1093826e-28;         // Electron mass
constexpr double MP          = 1.67262171e-24;        // Proton mass
constexpr double MN          = 1.67492728e-24;        // Neutron mass
constexpr double HPL         = 6.6260693e-27;         // Planck constant
constexpr double HBAR        = HPL / (2. * M_PI);     // Reduced Planck constant
constexpr double KBOL        = 1.3806505e-16;         // Boltzmann constant
constexpr double GNEWT       = 6.6742e-8;             // Gravitational constant
constexpr double SIG         = 5.670400e-5;           // Stefan-Boltzmann constant
constexpr double AR          = 4 * SIG / C_LIGHT_CGS; // Radiation constant
constexpr double THOMSON     = 0.665245873e-24;       // Thomson cross section
constexpr double COULOMB_LOG = 20.;                   // Coulomb logarithm
constexpr double ALPHAFS     = 0.007299270073;        // Fine structure constant ~ 1./137.
constexpr double GFERM       = 1.435850814e-49;       // Fermi constant
constexpr double GA          = -1.272323;             // Axial-vector coupling
constexpr double GA2         = GA * GA;               // Axial-vector coupling squared
constexpr double S2THW       = 0.222321;              // sin^2(Theta_W), Theta_W = Weinberg angle
constexpr double S4THW       = S2THW * S2THW;         // sin^4(Theta_W), Theta_W = Weinberg angle
constexpr double NUSIGMA0    = 1.7611737037e-44;      // Fundamental neutrino cross section
constexpr double AVO         = 6.0221417930e23;       // Avogadro's number
constexpr double AMU         = 1.66053878283e-24;     // Atomic Mass Unit

// Unit Conversion factors
constexpr double EV   = 1.60217653e-12;   // Electron-volt
constexpr double MEV  = 1.0e6 * EV;       // Mega-Electron-Volt
constexpr double GEV  = 1.0e9 * EV;       // Giga-Electron-Volt
constexpr double JY   = 1.e-23;           // Jansky
constexpr double PC   = 3.085678e18;      // Parsec
constexpr double AU   = 1.49597870691e13; // Astronomical unit
constexpr double YEAR = 31536000.;        // year in seconds
constexpr double DAY  = 86400.;           // day in seconds
constexpr double HOUR = 3600.;            // hour in seconds

// Macros
// ----------------------------------------------------------------------

// Primitive and conserved variables
constexpr int RHO = 0;
constexpr int UU = 1;
constexpr int U1 = 2;
constexpr int U2 = 3;
constexpr int U3 = 4;
constexpr int B1 = 5;
constexpr int B2 = 6;
constexpr int B3 = 7;
constexpr int NVAR_BASE = B3 + 1;

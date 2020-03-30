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
 * @file user.h
 * @author Julien Loiseau
 * @date April 2017
 * @brief User define for dimension and type 
 */

#ifndef _user_h_
#define _user_h_

#define OUTPUT
#define INTERNAL_ENERGY

static const size_t gdimension = EXT_GDIMENSION;
using type_t = double;

// FMM Taylor expansion order
#define fmm_order 3

// fix sph_kernel at compile time
// #define sph_kernel wendland_c4

// fix sph_viscosity at compile time
// #define sph_viscosity visc_constant

#endif // _user_h_

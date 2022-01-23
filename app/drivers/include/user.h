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

#define USER_H_STR_HELPER(x) #x
#define USER_H_STR(x) USER_H_STR_HELPER(x)
#define USER_H_DO_EXPAND(VAL)  VAL ## 1
#define USER_H_EXPAND(VAL)     USER_H_DO_EXPAND(VAL)

#define OUTPUT
#define INTERNAL_ENERGY

// Uncomment the next line to fix sph_kernel at compile time and
// enable vectorization:
// #define sph_kernel wendland_c4

static const size_t gdimension = EXT_GDIMENSION;
using type_t = double;

// FMM Taylor expansion order
#define fmm_order 3

// fix sph_kernel at compile time
// #define sph_kernel wendland_c4

// fix sph_viscosity at compile time
// #define sph_viscosity visc_constant

#if defined(EXT_EOS_TYPE) && (USER_H_EXPAND(EXT_EOS_TYPE) != 1)
// fix eos_type at compile time
#  define eos_type EXT_EOS_TYPE
#  pragma message "set at compile time: eos_type = " USER_H_STR(EXT_EOS_TYPE)
#endif

#endif // _user_h_

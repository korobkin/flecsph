/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2019 Triad National Security, LLC
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
 * @file utils.h
 * @authore Hyun Lim
 * @date Oct 2019
 * @brief Utilities for star tracker and GW radiation
 */

#ifndef _utils_h_
#define _utils_h_

#include "body.h"
#include "params.h"
#include "space_vector.h"
#include "tree.h"
#include "utils.h"
#include <cmath>

#define M_SUN_CGS 1.98847e33 // Solar mass in CGS
#define C_LIGHT_CGS 2.99792458e10 // Speed of light in CGS
constexpr double C_LIGHT = 3.424759e+2;
constexpr double C_LIGHT_NAT = 1.0;

#if 1

/*
 * A struct for tracking position
 * of star
 */
#define NDIMS 3 // number of dimensions
struct StarData {
  double mass;
  double mass1;
  double mass2;
  double radius1;
  double radius2;
  point_t center_of_mass;
  point_t ang_spin; // angular momentum vector
  point_t velocity;
};
typedef struct StarData StarData_t;

/*
 * A struct for tracking a binary system
 * with total COM and offsets
 */
#define NSTARS 2 // number of stars
struct BinaryData {
  point_t com;
  point_t velocity; // for COM
  point_t acc_gwcom[NSTARS]; // acceleration per star
  point_t offset[NSTARS]; // r_star - r_com
  double offset_norm[NSTARS]; // ||r_star - r_com||
  point_t ang_spin;
  // Backreaction for each star
  double separation;
  double total_mass;
  double reduced_mass; // dimensionaless reduced mass
  double omega_sq_correction; // omega^2 without Keplerian term
};
typedef struct BinaryData BinaryData_t;

#endif

#endif

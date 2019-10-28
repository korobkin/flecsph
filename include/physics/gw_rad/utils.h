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

#include <cmath>
#include "param.h"
#include "body.h"

#if 1

/*
 * A struct for tracking position
 * of star
 */
#define C_LIGHT 3.424759e+2
#define C_LIGHT_CGS 2.99792e10
#define C_LIGHT_NAT 1.0
#define NDIMS 3 //number of dimensions
typedef struct {
    double mass1;
    double mass2;
    double radius1;
    double radius2;
    double center_of_mass[NDIMS];
    double ang_spin[NDIMS]; //angular momentum vector
    double velocity[NDIMS];
} StarData_t;

/*
 * A struct for tracking a binary system
 * with total COM and offsets
 */
#define NSTARS 2 //number of stars
typedef struct {
   double center_of_mass[NDIMS];
   double velocity[NDIMS]; //for COM
   double acc_gwcom[NSTARS][NDIMS]; //acceleration per star
   double offset[NSTARS][NDIM]; //r_star - r_com
   double offset_norm[NSTARS]; // ||r_star - r_com|| 
   double ang_spin[NDIMS];
   // Backreaction for each star
   double separation;
   double total_mass;
   double reduced_mass;//dimensionaless reduced mass
   double omega_sq_correction; //omega^2 without Keplerian term
} BinaryData_t;

#endif

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
 * @file star_tracker.h
 * @authore Hyun Lim
 * @date Oct 2019
 * @brief Star tracking for PN 
 */

#if 1

#ifndef STAR_TRACKER_H
#define STAR_TRACKER_H

#include "utils.h"

// Get vector distance
double get_vec_dist(const double *v1, const double *v2) {
    double temp[NDIMS];
    double norm_diff;
    for (int i=0;i<NDIMS;i++){
      temp[i] = v1[i] - v2[i];
    }
    norm_diff = sqrt(temp[0]*temp[0]+temp[1]*temp[1]+temp[2]*temp[2]);
    return norm_diff;
{



#endif // STAR_TRACKER_H

#endif

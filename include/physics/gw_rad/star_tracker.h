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

/** 
 * Get vector distance
 * by the norm of two vectors ||v1 - v2||
 * Assume v1 and v2 are of length of NDIMS
 */
double get_vec_dist(const double *v1, const double *v2) {
    double temp[NDIMS];
    double norm_diff;
    for (int i=0;i<NDIMS;i++){
      temp[i] = v1[i] - v2[i];
    }
    norm_diff = sqrt(temp[0]*temp[0]+temp[1]*temp[1]+temp[2]*temp[2]);
    return norm_diff;
{

/** 
 * Find the maximum density for particles that lie within
 * radius of center of star
 */

void
find_max_density(std::vector<body>& bodies, double &maxrho1, double &maxrho2){

    maxrho1 = 0.0;
    maxrho2 = 0.0;
    for (auto b:bodies) {
        if(b->state() == STAR1) {
            maxrho1=std::max(maxrho1,b->getDensity());
        } else if (b->state() == STAR2) {
            maxrho2=std::max(maxrho2,b->getDensity());
        }
    }
    mpi_utils::reduce_max(maxrho1);
    mpi_utils::reduce_max(maxrho2);

}

void
cal_reduced_mass(std::vector<body>& bodies)
{

    double mass1, mass2
    for (auto b:bodies){
      if(b->state()==STAR1){
          mass1 += b->getmass()
      }
     }
      mpi_utils::reduce_sum(mass1);
      mpi_utils::reduce_sum(mass2);
      double reduced_mass = mass1*mass2/(mass1+mass2)

    

}

void
find_star_com(std::vector<body>& bodies, StarData_t* star, double &maxrho1, double &maxrho2){



}

void 
find_ang

void
omega_sq_kep(){

}

void meta_func(std::vector<body> &bodies, StarData_t* star, BinarySystem_t* bin)
{
    double maxrho1, maxrho2;
    find_max_density(bodies, maxrho1, maxrho2)
    find_com(bodies, star, maxrho1, maxrho2)
    find_ang
    cal_omega
    cal_kep
}

#endif // STAR_TRACKER_H

#endif

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
 * @file gw_rad.h
 * @authore Hyun Lim
 * @date Oct 2019
 * @brief Gravitatioanl radiation reaction via PN correction.
 * 	  Star tracking is done by adding "state" in body.
 * 	  Required information is calculated using the field
 */

#if 1

#ifndef GW_RAD_H
#define GW_RAD_H

#include "utils.h"

/* Star tracking part
 * First, we need to track the star by finding 
 * max density of the star and get COM of the system
 */


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
    std::cout<<"Maximum density for first star:"<<maxrho1<<std::endl;
    std::cout<<"Maximum density for second star:"<<maxrho2<<std::endl;
    mpi_utils::reduce_max(maxrho1);
    mpi_utils::reduce_max(maxrho2);

}

void
find_star_com(std::vector<body>& bodies, StarData_t* star){

  double global_max_density;
  double vector_distance;
  double total_mass = 0.0;
  double global_total_mass = 0.0;
  double ang_mom_star = 0.0;
  double ang_mom_part = 0.0;
  double momentum = 0.0;
  double global_momentum = 0.0;
  
  for (auto b:bodies){
    vector_distance = get_vec_dist(b->pos,max_density_part->pos);
    if(vector_distance < star->radius) 

    } 
  }


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

/*
 * Gravitational wave radiation-reaction part
 * We compute gw radition reactoin
 */

void
precompute_binary_system_props(BinaryData_t* system, 
			       const StarData_t* stars,
			       double param::gravitational_constant){
}

void
compute_particle_gw_acc(std::vector<body>& bodies,
                        const BinaryData_t* system,
	 	  	const StarData_t* stars,
 			const double* a_part_cart){
}

/*
 * Meta function that contains all functions from above.
 * This routine will be used to apply the calculations
 * to bodies.
 */




void gw_rad_PN(std::vector<body> &bodies, StarData_t* star, BinarySystem_t* bin)
{
    double maxrho1, maxrho2;
    find_max_density(bodies, maxrho1, maxrho2)
}

#endif // GW_RAD

#endif

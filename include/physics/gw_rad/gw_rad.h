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
find_max_density(std::vector<body>& bodies, body* mdp1, body* mdp2){

    double maxrho1 = 0.0;
    double maxrho2 = 0.0;
    for (auto b:bodies) {
        if (b->state() == STAR1) {
            if (maxrho1 < b->getDensity()){
	         maxrho1 = b->getDensity();
		 mdp1 = b;
	    }
        } else if (b->state() == STAR2) {
            if (maxrho2 < b->getDensity()){
	         maxrho2 = b->getDensity();
	 	 mdp2 = b;
            }
        }
    }
    std::cout<<"Maximum density for first star:"<<maxrho1<<std::endl;
    std::cout<<"Maximum density for second star:"<<maxrho2<<std::endl;
    mpi_utils::reduce_max(maxrho1);
    mpi_utils::reduce_max(maxrho2);

}

void
find_star_com(std::vector<body>& bodies, StarData_t* star){
  
  // Define quantities
  double vector_distance;

  body* mdp1;
  body* mdp2;

  double radius1 = 0.0, radius2 = 0.0;
  double momentum1=0.0, momentum2 = 0.0;
  double ang_mom_part1 = 0.0, ang_mom_part2 = 0.0;
  double ang_mom_star1 = 0.0, ang_mom_star2 = 0.0;
  
  point_t com1,com2;
  
  double total_mass1=0.0, total_mass2=0.0;
  
  //Find max density location. (We basically assume that this should be center of the star)
  find_max_density(bodies,mdp1,mdp2);

  for (auto b:bodies){
      if(b->state()==STAR1){
        vector_distance = get_vec_dist(b->coordinates(),mdp1->coordiantes());
	radius1 = std::max(std::distance(b->coordiates(),mdp1->cooridnates()),radius1);
        if(vector_distance < radius1) {
		com1 += (b->mass())*(b->coordinates());
		momentum1 += (b->mass())*(b->getVelocity());
		total_mass1 += b->mass();	
	}
      } else if(b->state()==STAR2){
        vector_distance = get_vec_dist(b->coordinates(),mdp2->coordiantes());
	radius2 = std::max(std::distance(b->coordiates(),mdp2->cooridnates()),radius2);
        if(vector_distance < star->radius2) {
		com2 += (b->mass())*(b->coordinates());
		momentum2 += (b->mass())*(b->getVelocity());
		total_mass2 += b->mass();
	}	
      }
    } 
     com1 /= total_mass1;
     com2 /= total_mass2;

  point_t global_momentum;
  double global_mass;

  for(size_t i = 0; i < bodies.size(), i++){
	global_momentum +=(bodies[i].mass())*(bodies[].getVelocity());
	global_mass +=bodies[i].mass()
  }
  mpi_utils::reduce_sum(global_momentum);
  mpi_utils::reduce_sum(global_mass);
  star->velocity = global_momentum;
  star->velocity /=global_mass;

  //Evaluate spin
   
  double particle_offset1, particle_offset2;  
  double rel_lin_mom1, rel_lin_mom2;  
  

  for (auto b:bodies) {
    if (b->state() == STAR1){
      particle_offset1 = b->coordinates() - com1;
      rel_lin_mom1 = b->velocity() - star->velocity;
      rel_lin_mom1 *= b->mass();
      ang_mom_part1 = vec_cross(particle_offset1, rel_lin_mom1);
      ang_mom_star1 += ang_mom_part1;
    else if (b->state() == STAR2){
      particle_offset2 = b->coordinates() - com2;
      rel_lin_mom2 = b->velocity() - star->velocity;
      rel_lin_mom2 *= b->mass();
      ang_mom_part1 = vec_cross(particle_offset2, rel_lin_mom1);
      ang_mom_star2 += ang_mom_part2;
    }      
    }
  }
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

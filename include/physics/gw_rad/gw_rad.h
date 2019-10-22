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

     mpi_utils::reduce_sum(com1);
     mpi_utils::reduce_sum(momentum1);
     mpi_utils::reduce_sum(total_mass1);

     mpi_utils::reduce_sum(com2);
     mpi_utils::reduce_sum(momentum2);
     mpi_utils::reduce_sum(total_mass2);

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
    mpi_utils::reduce_sum(ang_mom_star1);
    mpi_utils::reduce_sum(ang_mom_star2);
}


/*
 * Gravitational wave radiation-reaction part
 * We compute gw radition reactoin
 */

/*
 * Star velocity should be evaluated wrt to the velocity of
 * the COM
 */

void
precompute_binary_system_props(std::vector<body>& bodies,
			       BinaryData_t* system, 
			       const StarData_t* star){

  double star_mass1 = 0.0; star_mass2 = 0.0;

  for (auto b:bodies){
      if(b->state()==STAR1){
        if(vector_distance < radius1) {
		star_mass1 += b->mass();	
	}
      } else if(b->state()==STAR2){
        if(vector_distance < star->radius2) {
		star_mass2 += b->mass();
	}	
      }
    } 
  system->total_mass = star_mass1 + star_mass2;

  //Expansion parameter for PN theory (2.5PN order)
  double pn_param = param::gravitational_constant*(system->total_mass)
                    /(system->separation); //TODO : HL, need to check unit system
  //Omega term (related with orbital frequency) 
  //from generalization of the Kepler 3rd law
  double omega_sqd_kep = (system->total_mass)/(system->separation)
                         *(param::gravitational_constant/(system->separation))
                         /(system->separation);

  //Dimensionaless reduced mass
  system->reduced_mass = star_mass1/(system->total_mass)
		         *(star_mass2/(system->total_mass);
 
  //Binary system COM
  system->com = (star_mass1/(system->total_mass))*com1 
              + (star_mass2/(system->total_mass))*com2;

  //Relative velocities which are calculated wrt to the velocity of COM
  double rel_vel[NSTARS][NDIMS];

  system->velocity = star_mass1*star[0].velocity + star_mass2*star[1].velocity;
  system->velocity /= system->total_mass;

  double a_gwcm[NDIMS];
  // TODO : make parameter
  if (use_polar_coords){
    a_gwcm[0] = -1.*(system->omega_sq_correction/system->total_mass)*system->separation;
    a_gwcm[1] = -((32./5.)*pow(param::gravitational_constant,(7./2.))
		     *pow(system->total_mass,(5./2.))
                     *system->reduced_mass
                     /(pow(C_LIGHT_CGS,5.)*pow(system->separation,(9./2.))));
    a_gwcm[2] = 0.0;
    
    system->a_gwcm[0] = star_mass1*a_gwcm;
    system->a_gwcm[1] = star_mass2*a_gwcm;

    std::cout<<"Star COM acceleration"<<std::endl;
    std::cout<<"a_r_0 = "  <<a_gwcm[0][0]<<std::endl;
    std::cout<<"a_tan_0 = "<<a_gwcm[0][1]<<std::endl;
    std::cout<<"a_r_1 = "  <<a_gwcm[1][0]<<std::endl;
    std::cout<<"a_tan_1 = "<<a_gwcm[1][1]<<std::endl;

  } else if (use_vel_pos_basis) {
    double pre_factor = -(32./5.)*pow(pn_param,3.)
                        *(C_LIGHT_CGS/(system->separation))
                        *system->reduce_mass;
    std::cout<<"Velocity prefactor : "<<pre_factor<<std::endl;
    for (int i = 0; i < NSTARS; ++i) {
        system->a_gwcm[i] = - system->omega_sq_corrections*system->offset[i]
                            + pre_factor*rel_vel[i];

    }

  
  }
   
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

void 
gw_rad_PN(std::vector<body> &bodies, StarData_t* star, BinarySystem_t* system
               const double* a_part_cart)
{
    body *mdp1, *mdp2;
    find_max_density(bodies, mdp1, mdp2);
    fins_star_com(bodies, star);
    precompute_binary_system_props(bodies, system, star);
    compute_particle_gw_acc(bodies, system, star, a_part_cart);
}


if (enable_evaluate_gw_waveform) {
/*
 * Extract GW information via angle-averaged strain values
 * for + and x polarization
 * TODO : HL : I need to confirm this...
 */

 void 
 extract_gw_waveform(std::vector<body>) {

   double hp = 0.0;
   double hc = 0.0;

   //TODO : add quadrupole formula based on PN expansion
   // Ref : Zhuge et al. PRD.50.6247, 1994
   //       Blanchet. LRR-2014-2
   

 } //Evaluate GW waveform 

 /*
  * GW output
  * Ouputs for GW information that was calculated in previous
  */

 void
 gw_waveform_output(body_system<double,gdimension>& bs, const int rank) {
  
   static bool first_time = true;
   if (param::out_scalar_every <=0 ||
       physics::iteration % param::out_scalar_every !=0)
       return;

   // Compute GW information
   bs.get_all(extract_gw_radiation)l

   // output only from rank 0
   if (rank !=0) return;
   const char *filename = "PN_gw_info.dat";

   if (first_time) {
     // Generate output header
     /* HL :  Here, I assume that we only accept 3 dimensional case.
      *       It system of dimension is less than 3, code will be stopped
      * TODO : will be generalized once we have lower dimensional case
      *         such as axisymmetry case
      */
     switch(gdimension){
     case 1:
       std::cout<<"System of dimension must be 3"<<std::endl;
       assert(false);
     break;

     case 2:
       std::cout<<"System of dimension must be 3"<<std::endl;
       assert(false);
     break;
 
     case 3:
     default:
       oss_header
         << "# GW Waveform Data:"<<std::endl;
         << "# 1:iteration 2:time 3:timestep"<<std::endl
         << "# 4:rh+ 5:rhx">>std::endl;
     }
    
     std::ofstream out(filename);
     out << oss_header.str();
     out.close();
     first_time = false;
   }
 
   std::ostringstream oss_data;
   oss_data << std::setw(14) << physics::iteration
      << std::setw(24) << std::scientific << std::setprecision(22)
      << physics::totaltime << std::setw(20) << physics::dt << " "
      << hp << " " << hc << " "<<std::endl;
  
   // Open file in append mode
   std::ofstream out(filename,std::ios_base::app);
   out << oss_data.str();
   out.close;

 }// gw_waveform_output

}// enable_evaluate_gw_waveform

#endif // GW_RAD

#endif

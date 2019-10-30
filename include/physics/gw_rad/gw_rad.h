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
    norm_diff = std::sqrt(temp[0]*temp[0]+temp[1]*temp[1]+temp[2]*temp[2]);
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

  for (int i = 0; i < NSTARS; ++i){
      rel_vel[i] = vec_diff(star[i].velocity, system->velocity);
  }

  for (int i = 0; i < NSTARS; ++i){
      // Star offsets from total com
      system->offset[i] = vec_diff(star[i].center_of_mass, system->center_of_mass);
      //Radius from common COM to COM of each star
      system->offset_norm[i] = std::sqrt(vec_dot(system->offset[i],system->offset[i]));
  }

  // Set spin
  for (int i = 0; i < NDIMS; ++i){
     system->ang_spin[i] = 0.0;
  }

  double star_L[NDIMS], star_P[NDIMS];
  
  for (int i = 0; i < NSTARS; ++i){
     star_P = star[i].mass*rel_vel[i];   
     star_L = vec_cross(system->offset[i],star_P);
     system->ang_spin += star_L;
  }

  //Some prefactor
  double red_mass_pre_fac = 6.+(41./4.)*(system->reduced_mass)
                            + (system->reduced_mass)*(system->reduced_mass);

  // Omega correction term
  system->omega_sq_correction = omega_sqd_kep*pn_param
                                * ((system->reduced_mass - 3.) 
                                + pn_param*red_mass_pre_fac);

  
  // Acceleration for COM in polar coordinates 
  // HL : Need to check
  if (use_polar_coords){
    double a_gwcm[NDIMS];
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
  
  } else {
    std::cout<<"Wrong choice"<<std::endl;
    assert(false);
  }
   
}

// Collection of functions that are useful
// TODO : we can clean this up

double get_angle(const double* offsets){
  return atan(offsets[1],offsets[0]);
}

double get_ang_vel(const double* offsets, const double* vels){
  double omega;
  omega = offsets[0]*vels[1] - offsets[1]*vels[0];
  omega /= (offsets[0]*offsets[0] + offsets[2]*offsets[2]);
  return omega;
}

template <typename T>
T get_sign(T qt){
  return (T(0) < val) - (T(0) > val);
}

// Now compute particle GW acceleration

void
compute_particle_gw_acc(std::vector<body>& bodies,
                        const BinaryData_t* system,
	 	  	const StarData_t* stars,
 			const double* a_part_cart){

  double particle_offset[NDIMS];
  double a_par_polar[NDIMS];

  double r, r_spherical;
  // Radius measured at the COM of my star
  double rcm
  // Angular info of the star
  double star_theta[NSTARS], star_omega[NSTARS];
  double rel_vel[NSTARS][NDIMS];

  // Particel offset : r_part - r_system_com
  particle_offset = b->coordinates() - system->center_of_mass;
  r_spherical = std::sqrt(vec_dot(particle_offset,particle_offset));

  if (use_polar_coords){
    for (int i = 0; i < NSTARS; ++i) {
      rel_vel[i] = vec_diff(stars[i].velocity, system->velocity);
      star_theta[i] = get_angle(system->offset[i]);
      star_omega[i] = get_ang_vel(system->offset[i],rel_vel[i]);
    }

    // Check the sign to have same rotation direction for both stars
    double sign_omega = get_sign(star_omega[0]);
    
    // Get distance from COM for each star
    if (b->state() == STAR1){
        rcm = system->offset_norm[0];
        r = vec_dot(particle_offset,system->offset[0]);
        r /= system->offset_norm[0];
    } else if (b->state() == STAR2){
        rcm = system->offset_norm[1];
        r = vec_dot(particle_offset,system->offset[1]);
        r /= system->offset_norm[1];
    } else {
       std::cout<<"This isn't the case. Set position to origin"<<std::endl;
       rcm = 0.0; r = 0.0;
    }

    if (b->state() == STAR1 || b->state() == STAR2){
        if (polar_radial_dependence){
           a_par_polar[0] = system->a_gwcm[b->state()][0];
           a_par_polar[1] = sign_omega*(r/rcm)
                         *system->a_gwcm[b->state()][1];
           a_par_polar[2] = 0.0;
        } else {
           // no radial dependence : unlocked rigid rotation
           a_par_polar[0] = system->a_gwcm[b->state()][0];
           a_par_polar[1] = sign_omega*system->a_gwcm[b->state()][1];
           a_par_polar[2] = 0.0;
        }

      // Acceleration for particle in Cartesian coordinate
      // We pass the acceleration values in polar coordinates
      // to Cartesian coordinates via the rotation endomorphism
      a_part_cart[0] = cos(star_theta[b->state()])*a_par_polar[0]
                      -sin(star_theta[b->state()])*a_par_polar[1];
      a_part_cart[1] = sin(star_theta[b->state()])*a_par_polar[0]
                      +cos(star_theta[b->state()])*a_par_polar[1];
      a_part_cart[2] = 0.0;
    } else {
      a_part_cart[0] = 0.0;
      a_part_cart[1] = 0.0;
      a_part_cart[2] = 0.0;
    }
 } else if (use_vel_pos_basis) {
 // Using velocity/position basis
    if (b->state() == STAR1 || b->state() == STAR2){
       a_part_cart[0] = system->a_gwcm[b->state()][0];
       a_part_cart[1] = system->a_gwcm[b->state()][1];
       a_part_cart[2] = system->a_gwcm[b->state()][2];
    } else {
       a_part_cart[0] = 0.0;
       a_part_cart[1] = 0.0;
       a_part_cart[2] = 0.0;
    }
 } else {
  std::cout<<"Wrong choice"<<std::endl:
  assert(false);
 }

} // compute_gw_particle_acc

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

/*
 * Extract GW information via angle-averaged strain values
 * for + and x polarization
 * TODO : HL : I need to confirm this...
 */

 void 
 extract_gw_waveform(body& particle, std::vector<body> &bodies) {

   //Define angle averaged value of strain:
   // <rh_+> and <rh_x>. 

   double strain_hp = 0.0; // Plus polarization
   double strain_hc = 0.0; // Cross polarization

   if (enable_evaluate_gw_waveform){
   // Ref : Zhuge et al. PRD.50.6247, 1994
   //       Blanchet. LRR-2014-2
   //       van den Broek et al. MNRAS 425, L24-L27, 2012
   // HL : Since we don't have tensor contribution, I will list 
   //      non_zero component of quadrupole moment.

     //TODO  : Observer position?
     // HL : We may not need this for below formulation so I commented out
     #if 0
     point_t obs_pos={1,1,1};

     for (auto b:bodies) {
      double r_dist = distance(b->coordinates(),obs_pos);
     
      double qxxc = b->coordinates()[0] - obs_pos[0];
      double qxyc = b->coordinates()[0] - obs_pos[1];
      double qxzc = b->coordinates()[1] - obs_pos[2];
      double qyyc = b->coordinates()[1] - obs_pos[1];
      double qyzc = b->coordinates()[1] - obs_pos[2];
      double qzzc = b->coordinates()[2] - obs_pos[2];

      double qxx = 3*.b->getMass()*qxxc*qxxc - r_dist;
      double qxy = 3*.b->getMass()*qxyc*qxyc;
      double qxz = 3*.b->getMass()*qxzc*qxzc;
      double qyy = 3*.b->getMass()*qyyc*qyyc - r_dist;
      double qyz = 3*.b->getMass()*qyzc*qyzc;
      double qzz = 3*.b->getMass()*qzzc*qzzc - r_dist;

     }
     
     mpi_utils::reduce_sum(qxx);
     mpi_utils::reduce_sum(qxy);
     mpi_utils::reduce_sum(qxz);
     mpi_utils::reduce_sum(qyy);
     mpi_utils::reduce_sum(qyz);
     mpi_utils::reduce_sum(qzz);

      // Express quadrupole moments in terms of orthonormal spherical coordinates
      
      // Get polar and azimuthal angle
      // HL : Here, I set both angles to be zero i.e. observer
      //      located on the axix.
      // TODO : Generalized this
      double polar_ang = 0.0;
      double az_ang = 0.0;
      
      double Ithetatheta = (qxx*cos(az_ang)*cos(az_ang) + qyy*sin(az_ang)*sin(az_ang)
                         + qxy*sin(2.0*az_ang))*cos(polar_ang)*cos(polar_ang)
                         + qzz*sin(polar_ang)*sin(polar_ang)
                         - (qxz*cos(az_ang)+qyz*sin(az_ang))*sin(2.0*polar_ang);
      double Iphiphi = qzz*sin(az_ang)*sin(az_ang) + qyy*cos(az_ang)*cos(az_ang)
                     - qxy*sin(2.0*az_ang);
      double Ithetaphi = 0.5*(qyy-qxx)*cos(polar_ang)*sin(2.0*az_ang) 
                       + qxy*cos(polar_ang)*cos(2.0*az_ang)
                       + (qxz*sin(az_ang) - qyz*cos(az_ang))*sin(polar_ang);
     
     // Time derivatives for moments to compute strain
     // HL : How?
     double Ithetatheta_dtdt = 0.0;
     double Iphiphi_dtdt = 0.0;
     double Ithetaphi_dtdt = 0.0;


      // Now compute strain for waveforms
      // TODO : units?
      strain_hp = param::gravitational_constant/(pow(C_LIGHT_CGS,4.0)*r_dist)
                  *(Ithetatheta_dtdt - Iphiphi_dtdt);
      strain_hc = (2.0*param::gravitational_constant)/(pow(C_LIGHT_CGS,4.0)*r_dist)
                  *(Ithetaphi_dtdt);
               
      #endif
      // HL : Here, I propose some alternative way to compute rhx and rh+
      //      This is based on just using particles' position, velocity, 
      //      and acceleration to compute quadrupole moments and its derivatives

      for (auto b:bodies){
        // Compute second time derivatives of each components of quadrupole moments

        double qxx_dtdt = 0.0, qyy_dtdt = 0.0; qzz_dtdt = 0.0;
        double qxy_dtdt = 0.0, qxz_dtdt = 0.0; qyz_dtdt = 0.0;

        qxx_dtdt = 2.0/3.0 * b->getMass() 
                 * (2.0*b->coordinates()[0]*b->getAcceleration()[0]
                    - b->coordinates()[1]*b->getAcceleration()[1] 
                    - b->coordinates()[2]*b->getAcceleration()[2]
                    + 2.0*b->getVelocity()[0]*b->getVelocity()[0]
                    - b->getVelocity()[1]*b->getVelocity()[1]
                    - b->getVelocity()[2]*b->getVelocity()[2]
                   );

        qxy_dtdt = b->getMass()
                 * (b->coordinats()[0]*b->getAcceleration()[1] 
                    + b->coordinates()[1]*b->getAcceleration()[0]
                    + 2.0*b->getVelocity()[0]*b->getVelocity()[1]
                   );

        qxz_dtdt = b->getMass()
                 * (b->coordinats()[0]*b->getAcceleration()[2] 
                    + b->coordinates()[2]*b->getAcceleration()[0]
                    + 2.0*b->getVelocity()[0]*b->getVelocity()[2]
                   );

        qyy_dtdt = 2.0/3.0 * b->getMass() 
                 * (2.0*b->coordinates()[1]*b->getAcceleration()[1]
                    - b->coordinates()[0]*b->getAcceleration()[0] 
                    - b->coordinates()[2]*b->getAcceleration()[2]
                    + 2.0*b->getVelocity()[1]*b->getVelocity()[1]
                    - b->getVelocity()[0]*b->getVelocity()[0]
                    - b->getVelocity()[2]*b->getVelocity()[2]
                   );
        
        qyz_dtdt = b->getMass()
                 * (b->coordinats()[1]*b->getAcceleration()[2] 
                    + b->coordinates()[2]*b->getAcceleration()[1]
                    + 2.0*b->getVelocity()[1]*b->getVelocity()[2]
                   );

        qzz_dtdt = 2.0/3.0 * b->getMass() 
                 * (2.0*b->coordinates()[2]*b->getAcceleration()[2]
                    - b->coordinates()[1]*b->getAcceleration()[1] 
                    - b->coordinates()[0]*b->getAcceleration()[0]
                    + 2.0*b->getVelocity()[2]*b->getVelocity()[2]
                    - b->getVelocity()[1]*b->getVelocity()[1]
                    - b->getVelocity()[0]*b->getVelocity()[0]
                   );

        // Using symmetric property defining remaining values
        double qyx_dtdt = qxy_dtdt;
        double qzx_dtdt = qxz_dtdt;
        double qzy_dtdt = qyz_dtdt;
      }

     mpi_utils::reduce_sum(qxx_dtdt);
     mpi_utils::reduce_sum(qxy_dtdt);
     mpi_utils::reduce_sum(qxz_dtdt);
     mpi_utils::reduce_sum(qyy_dtdt);
     mpi_utils::reduce_sum(qyz_dtdt);
     mpi_utils::reduce_sum(qzz_dtdt);

     // Compute angle averaged strain values
     double strain_hp_sq = 0.0, strain_hc_sq = 0.0;

     strain_hp_sq = 4.0/15.0*((qxx_dtdt - qzz_dtdt)*(qxx_dtdt - qzz_dtdt)
                              + (qyy_dtdt - qzz_dtdt)*(qyy_dtdt - qzz_dtdt)
                              + qxz_dtdt*qxz_dtdt + qyz_dtdt*qyz_dtdt)
                  + 1.0/10.0*(qxx_dtdt - qyy_dtdt)*(qxx_dtdt - qyy_dtdt)
                  + 14.0/15.0*qxy_dtdt*qxy_dtdt;
     
     strain_hc_sq = 1.0/6.0*(qxx_dtdt - qyy_dtdt)*(qxx_dtdt - qyy_dtdt)
                  + 2.0/3.0*qxy_dtdt*qxy_dtdt 
                  + 4.0/3.0*(qxz_dtdt*qxz_dtdt + qyz_dtdt*qyz_dtdt);

     strain_hp = std::sqrt(starin_hp_sq);
     strain_hc = std::sqrt(strain_hc_sq);

     #if 0
     // HL : I just added this for future ref
     if(do_2p5PN){
     //TODO : add higher order PN
     }
     #endif

   } 
   // If we don't include this routine, we will
   // only have zeros for strain TODO : Good? 

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
   bs.get_all(extract_gw_radiation);

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
      << strain_hp << " " << starin_hc << " "<<std::endl;
  
   // Open file in append mode
   std::ofstream out(filename,std::ios_base::app);
   out << oss_data.str();
   out.close;

 }// gw_waveform_output

#endif // GW_RAD

#endif

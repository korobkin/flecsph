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
 *        Star tracking is done by adding "state" in body.
 *        Required information is calculated using the field
 */

#ifndef _GW_RAD_H_
#define _GW_RAD_H_

#include "utils_gw.h"

/* Star tracking part
 * First, we need to track the star by finding
 * max density of the star and get COM of the system
 */

/**
 * Find the maximum density for particles that lie within
 * radius of center of star
 */

void
find_max_density(std::vector<body>& bodies, body* mdp1, body* mdp2){

    double maxrho1 = 0.0;
    double maxrho2 = 0.0;
    for (auto b:bodies) {
        if (b.state() == STAR1) {
            if (maxrho1 < b.getDensity()){
             maxrho1 = b.getDensity();
         mdp1 = &b;
        }
        } else if (b.state() == STAR2) {
            if (maxrho2 < b.getDensity()){
             maxrho2 = b.getDensity();
          mdp2 = &b;
            }
        }
    }
  std::cout << "Maximum density for first star:" << maxrho1 << std::endl;
  std::cout << "Maximum density for second star:" << maxrho2 << std::endl;
  mpi_utils::reduce_max(maxrho1);
  mpi_utils::reduce_max(maxrho2);
}

void
find_star_com(std::vector<body>& bodies, StarData_t* star){

  // Define quantities
  double vector_distance;

  body *mdp1{NULL}, *mdp2{NULL};

  double radius1 = 0.0, radius2 = 0.0;
  point_t momentum1 = 0.0;
  point_t momentum2 = 0.0;
  point_t ang_mom_part1 = 0.0, ang_mom_part2 = 0.0;
  point_t ang_mom_star1 = 0.0, ang_mom_star2 = 0.0;

  point_t com1, com2;

  double total_mass1 = 0.0, total_mass2 = 0.0;

  // Find max density location. (We basically assume that this should be center
  // of the star)
  find_max_density(bodies, mdp1, mdp2);

  for (auto b:bodies){
    if (b.state()==STAR1){
      vector_distance = distance(b.coordinates(), mdp1->coordinates());
      radius1 = std::max(distance(b.coordinates(),mdp1->coordinates()), radius1);
      if(vector_distance < radius1) {
        com1 += (b.mass())*(b.coordinates());
        momentum1 += (b.mass())*(b.getVelocity());
        total_mass1 += b.mass();
      }
    }
    else if (b.state()==STAR2) {
      vector_distance = distance(b.coordinates(), mdp2->coordinates());
      radius2 = std::max(distance(b.coordinates(),mdp2->coordinates()), radius2);
      if(vector_distance < star->radius2) {
        com2 += (b.mass())*(b.coordinates());
        momentum2 += (b.mass())*(b.getVelocity());
        total_mass2 += b.mass();
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

  for(auto b:bodies){
    global_momentum +=(b.mass())*(b.getVelocity());
    global_mass += b.mass();
  }
  mpi_utils::reduce_sum(global_momentum);
  mpi_utils::reduce_sum(global_mass);
  star->velocity = global_momentum;
  star->velocity /= global_mass;

  // Evaluate spin

  point_t particle_offset1, particle_offset2;
  point_t rel_lin_mom1, rel_lin_mom2;


  for (auto b:bodies) {
    if (b.state() == STAR1){
      particle_offset1 = b.coordinates() - com1;
      rel_lin_mom1 = b.getVelocity() - star->velocity;
      rel_lin_mom1 *= b.mass();
      ang_mom_part1 = cross(particle_offset1, rel_lin_mom1);
      ang_mom_star1 += ang_mom_part1;
    } else if (b.state() == STAR2){
      particle_offset2 = b.coordinates() - com2;
      rel_lin_mom2 = b.getVelocity() - star->velocity;
      rel_lin_mom2 *= b.mass();
      ang_mom_part1 = cross(particle_offset2, rel_lin_mom1);
      ang_mom_star2 += ang_mom_part2;
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
                   StarData_t* star){

  double star_mass1 = 0.0, star_mass2 = 0.0;

  double vector_distance;
  double radius1 = 0.0, radius2 = 0.0;

  body *mdp1{NULL}, *mdp2{NULL};

  point_t com1, com2;

  find_max_density(bodies, mdp1, mdp2);
  find_star_com(bodies, star);

  for (auto b:bodies){
    if(b.state()==STAR1){
      vector_distance = distance(b.coordinates(), mdp1->coordinates());
      radius1 = std::max(distance(b.coordinates(),mdp1->coordinates()),radius1);
      if(vector_distance < radius1)
        star_mass1 += b.mass();
    } 
    else if(b.state()==STAR2){
      vector_distance = distance(b.coordinates(), mdp2->coordinates());
      radius2 = std::max(distance(b.coordinates(),mdp2->coordinates()),radius2);
      if(vector_distance < radius2)
        star_mass2 += b.mass();
    }
  }
  system->total_mass = star_mass1 + star_mass2;

  // Expansion parameter for PN theory (2.5PN order)
  double pn_param =
    param::gravitational_constant * (system->total_mass) /
    (system->separation); // TODO : HL, need to check unit system
  // Omega term (related with orbital frequency)
  // from generalization of the Kepler 3rd law
  double omega_sqd_kep =
    (system->total_mass) / (system->separation) *
    (param::gravitational_constant / (system->separation)) /
    (system->separation);

  // Dimensionaless reduced mass
  system->reduced_mass =
    star_mass1 / (system->total_mass) * (star_mass2 / (system->total_mass));

  // Binary system COM
  system->com = (star_mass1 / (system->total_mass)) * com1 +
                (star_mass2 / (system->total_mass)) * com2;

  // Relative velocities which are calculated wrt to the velocity of COM
  point_t rel_vel[NSTARS];

  system->velocity =
    star_mass1 * star[0].velocity + star_mass2 * star[1].velocity;
  system->velocity /= system->total_mass;

  for(int i = 0; i < NSTARS; ++i) {
    rel_vel[i] = star[i].velocity - system->velocity;
  }

  for(int i = 0; i < NSTARS; ++i) {
    // Star offsets from total com
    system->offset[i] = star[i].center_of_mass - system->com;
    // Radius from common COM to COM of each star
    system->offset_norm[i] = magnitude(system->offset[i]);
  }

  // Set spin
  for(int i = 0; i < NDIMS; ++i) {
    system->ang_spin[i] = 0.0;
  }

  point_t star_L, star_P;

  for(int i = 0; i < NSTARS; ++i) {
    star_P = star[i].mass * rel_vel[i];
    star_L = cross(system->offset[i], star_P);
    system->ang_spin += star_L;
  }

  // Some prefactor
  double red_mass_pre_fac = 6. + (41. / 4.) * (system->reduced_mass) +
                            (system->reduced_mass) * (system->reduced_mass);

  // Omega correction term
  system->omega_sq_correction =
    omega_sqd_kep * pn_param *
    ((system->reduced_mass - 3.) + pn_param * red_mass_pre_fac);

  // Acceleration for COM in polar coordinates
  // HL : Need to check
  if(param::use_polar_coords) {
    point_t a_gwcm;
    a_gwcm[0] = -1. * (system->omega_sq_correction / system->total_mass) *
                system->separation;
    a_gwcm[1] = -((32. / 5.) * pow(param::gravitational_constant, (7. / 2.)) *
                  pow(system->total_mass, (5. / 2.)) * system->reduced_mass /
                  (pow(C_LIGHT_CGS, 5.) * pow(system->separation, (9. / 2.))));
    a_gwcm[2] = 0.0;

    system->acc_gwcom[0] = star_mass1 * a_gwcm;
    system->acc_gwcom[1] = star_mass2 * a_gwcm;

#if 0
    std::cout<<"Star COM acceleration"<<std::endl;
    std::cout<<"a_r_0 = "  <<a_gwcm[0][0]<<std::endl;
    std::cout<<"a_tan_0 = "<<a_gwcm[0][1]<<std::endl;
    std::cout<<"a_r_1 = "  <<a_gwcm[1][0]<<std::endl;
    std::cout<<"a_tan_1 = "<<a_gwcm[1][1]<<std::endl;
#endif
  }
  else if(param::use_vel_pos_basis) {
    double pre_factor = -(32. / 5.) * pow(pn_param, 3.) *
                        (C_LIGHT_CGS / (system->separation)) *
                        system->reduced_mass;
    std::cout << "Velocity prefactor : " << pre_factor << std::endl;
    for(int i = 0; i < NSTARS; ++i) {
      system->acc_gwcom[i] = -system->omega_sq_correction * system->offset[i] +
                             pre_factor * rel_vel[i];
    }
  }
  else {
    std::cerr << "Wrong choice" << std::endl;
    assert(false);
  }
}

// Collection of functions that are useful
// TODO : we can clean this up

double
get_angle(point_t offsets) {
  return atan2(offsets[1], offsets[0]);
}

double
get_ang_vel(point_t offsets, point_t vels) {
  double omega;
  omega = offsets[0] * vels[1] - offsets[1] * vels[0];
  omega /= (offsets[0] * offsets[0] + offsets[2] * offsets[2]);
  return omega;
}

template<typename T>
T
get_sign(T val) {
  return (T(0) < val) - (T(0) > val);
}

// Now compute particle GW acceleration

void
compute_particle_gw_acc(std::vector<body>& bodies,
    BinaryData_t* system,
    StarData_t* stars,
    double* a_part_cart)
{

  point_t particle_offset;
  point_t a_par_polar;

  double r, r_spherical;
  // Radius measured at the COM of my star
  double rcm;
  // Angular info of the star
  double star_theta[NSTARS], star_omega[NSTARS];
  point_t rel_vel[NSTARS];

  // Particel offset : r_part - r_system_com
  for(auto b:bodies){
    particle_offset = b.coordinates() - system->com;
  }
  r_spherical = magnitude(particle_offset);

  if(param::use_polar_coords) {
    for(int i = 0; i < NSTARS; ++i) {
      rel_vel[i] = stars[i].velocity - system->velocity;
      star_theta[i] = get_angle(system->offset[i]);
      star_omega[i] = get_ang_vel(system->offset[i], rel_vel[i]);
    }

    // Check the sign to have same rotation direction for both stars
    double sign_omega = get_sign(star_omega[0]);

    // Get distance from COM for each star
    for (auto b:bodies){
      if (b.state() == STAR1){
        rcm = system->offset_norm[0];
        r = dot(particle_offset, system->offset[0]);
        r /= system->offset_norm[0];
      } 
      else if (b.state() == STAR2){
        rcm = system->offset_norm[1];
        r = dot(particle_offset, system->offset[1]);
        r /= system->offset_norm[1];
      }
      else {
        std::cout << "This isn't the case. Set position to origin" << std::endl;
        rcm = 0.0;
        r = 0.0;
      }
    }

    for (auto b:bodies){
      if (b.state() == STAR1 || b.state() == STAR2){
        auto st = (b.state() == STAR1) ? 0 : 1; 
        if (param::polar_radial_dependence){
           a_par_polar[0] = system->acc_gwcom[st][0];
           a_par_polar[1] = sign_omega*(r/rcm)
                         *system->acc_gwcom[st][1];
           a_par_polar[2] = 0.0;
        } 
        else {
           // no radial dependence : unlocked rigid rotation
           a_par_polar[0] = system->acc_gwcom[st][0];
           a_par_polar[1] = sign_omega*system->acc_gwcom[st][1];
           a_par_polar[2] = 0.0;
        }

        // Acceleration for particle in Cartesian coordinate
        // We pass the acceleration values in polar coordinates
        // to Cartesian coordinates via the rotation endomorphism
        a_part_cart[0] = cos(star_theta[st])*a_par_polar[0]
                    -sin(star_theta[st])*a_par_polar[1];
        a_part_cart[1] = sin(star_theta[st])*a_par_polar[0]
                    +cos(star_theta[st])*a_par_polar[1];
        a_part_cart[2] = 0.0;
      }
      else {
        a_part_cart[0] = 0.0;
        a_part_cart[1] = 0.0;
        a_part_cart[2] = 0.0;
      }
   }
 } 
 else if (param::use_vel_pos_basis) {
 // Using velocity/position basis
   for (auto b:bodies){
     if (b.state() == STAR1 || b.state() == STAR2){
       auto st = (b.state() == STAR1) ? 0 : 1; 
       a_part_cart[0] = system->acc_gwcom[st][0];
       a_part_cart[1] = system->acc_gwcom[st][1];
       a_part_cart[2] = system->acc_gwcom[st][2];
     } else {
       a_part_cart[0] = 0.0;
       a_part_cart[1] = 0.0;
       a_part_cart[2] = 0.0;
     }
   }
 } else {
   std::cerr<<"Wrong choice"<<std::endl;
   assert(false);
 }

} // compute_gw_particle_acc

/*
 * Meta function that contains all functions from above.
 * This routine will be used to apply the calculations
 * to bodies.
 */

void
gw_rad_PN(std::vector<body> &bodies)
{
  StarData_t star;
  BinaryData_t system;
  double a_part_cart;
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  body *mdp1{NULL}, *mdp2{NULL};
  find_max_density(bodies, mdp1, mdp2);
  find_star_com(bodies, &star);
  precompute_binary_system_props(bodies, &system, &star);
  compute_particle_gw_acc(bodies, &system, &star, &a_part_cart);
}
#endif //_GW_RAD_H_

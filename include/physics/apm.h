/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2020 Triad National Security, LLC
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
 * @file apm.h
 * @date December 2020
 * @brief Interface to select various artificial pressure computation
 */

#pragma once

#include "params.h"
#include "tree.h"
#include "user.h"
#include <boost/algorithm/string.hpp>
#include <math.h>
#include <stdlib.h>
#include "density_profiles.h"
#include "kernels.h"
#include "eos.h"

#define SQ(x) ((x) * (x))
#define CU(x) ((x) * (x) * (x))

//TODO : make it as parameter
//#define LE_density
//#define APM_MONITOR
//#define POS_CORRECTION

// TODO : check position corrector for all the other cases. HL : Need to make separate routine?

namespace apm {

template<param::apm_type_keyword K>
point_t compute_apm_acc(body & particle, std::vector<body *> &nbs);

typedef point_t (*compute_apm_acc_t)(body & particle, std::vector<body *> &nbs);

/**
 * @brief  Zero APM force i.e. do nothing
 */
template<>
point_t
compute_apm_acc<param::zero_apm>(body & particle, std::vector<body *> &nbs) {

    point_t a_apm = 0.0;

    return a_apm;

}
/**
 * @brief  kilonova spherical-ejecta 
 */
template<>
point_t
compute_apm_acc<param::kn_ejecta>(body & particle, std::vector<body *> &nbs) {
   using namespace param;

    // Relative density error
    double Pi_a = 0.0, Pi_b = 0.0;
    // Base pressure
    double P0 = apm_base_pressure;

    // Acceleration from artificial pressure
    point_t a_apm = 0.0;

    // this particle (index 'a')
    const double h_a = particle.radius(),
                 rho_a = particle.getDensity();
    const point_t pos_a = particle.coordinates();

    // neighbor particles (index 'b')
    const int n_nb = nbs.size();
    double rho_[n_nb], m_[n_nb], h_[n_nb];
    point_t pos_[n_nb], DiWa_[n_nb];
    
    // compute kernel gradient
    for (int b = 0; b < n_nb; ++b) {
      const body * const nb = nbs[b];
      rho_[b] = nb->getDensity();
      pos_[b] = nb->coordinates();
      h_[b] = nb->radius();
      m_[b] = nb->mass()*(pos_[b] != pos_a);
      const point_t pos_ab = pos_a - pos_[b];
      const double h_ab = 0.5*(h_a + h_[b]);
      DiWa_[b] = kernels::sph_kernel_gradient(pos_ab,h_ab);
    }

    //Kilonova profile
    const double r_a = flecsi::magnitude(pos_a);
    double rho_a_target = kn_ejecta_mass/CU(sphere_radius)*density_profiles::rho_kn_ejecta(r_a/sphere_radius);
           rho_a_target *= M_SUN_CGS;
    
    // compute apm acceleration and position corrector
    point_t dr_apm; // position updates
    for (int b = 0; b < n_nb; ++b) {
      const double r_b = flecsi::magnitude(pos_[b]);
      double rho_b_target = kn_ejecta_mass/CU(sphere_radius)*density_profiles::rho_kn_ejecta(r_b/sphere_radius);
             rho_b_target *= M_SUN_CGS;
      Pi_a = std::max(1.0 + (rho_a - rho_a_target)/rho_a_target, 0.1);
      Pi_b = std::max(1.0 + (rho_[b] - rho_b_target)/rho_b_target, 0.1);
      a_apm = -P0*m_[b] * (Pi_a + Pi_b)/(rho_a*rho_[b]) * DiWa_[b];
#ifdef POS_CORRECTION
      dr_apm = - apm_pos_prefactor*h_a*h_a*m_[b]*(Pi_a + Pi_b)/rho_[b] * DiWa_[b];
#endif
    }
    
    #ifdef POS_CORRECTION
    particle.set_coordinates(pos_a+dr_apm);
    #endif

    return a_apm;

}

/**
 * @brief  Sharp spherical density profile
 */

 //TODO : check this profile more
template<>
point_t
compute_apm_acc<param::sharp_spherical>(body & particle, std::vector<body *> &nbs) {
   using namespace param;

    // Relative density error
    double Pi_a = 0.0, Pi_b = 0.0;
    // Base pressure
    double P0 = apm_base_pressure;

    // Acceleration from artificial pressure
    point_t a_apm = 0.0;

    // this particle (index 'a')
    const double h_a = particle.radius(),
                 rho_a = particle.getDensity();
    const point_t pos_a = particle.coordinates();

    // neighbor particles (index 'b')
    const int n_nb = nbs.size();
    double rho_[n_nb], m_[n_nb], h_[n_nb];
    point_t pos_[n_nb], DiWa_[n_nb];
    
    // compute kernel gradient
    for (int b = 0; b < n_nb; ++b) {
      const body * const nb = nbs[b];
      rho_[b] = nb->getDensity();
      pos_[b] = nb->coordinates();
      h_[b] = nb->radius();
      m_[b] = nb->mass()*(pos_[b] != pos_a);
      const point_t pos_ab = pos_a - pos_[b];
      const double h_ab = 0.5*(h_a + h_[b]);
      DiWa_[b] = kernels::sph_kernel_gradient(pos_ab,h_ab);
    }

    const double r_a = flecsi::magnitude(pos_a);
    double rho_a_target = density_profiles::rho_sharp_spherical(r_a);
    
    // compute apm acceleration
    for (int b = 0; b < n_nb; ++b) {
      const double r_b = flecsi::magnitude(pos_[b]);
      double rho_b_target = density_profiles::rho_sharp_spherical(r_b);
      Pi_a = std::max(1.0 + (rho_a - rho_a_target)/rho_a_target, 0.1);
      Pi_b = std::max(1.0 + (rho_[b] - rho_b_target)/rho_b_target, 0.1);
      a_apm = -P0*m_[b] * (Pi_a + Pi_b)/(rho_a*rho_[b]) * DiWa_[b];
    }

    return a_apm;

}
/**
 * @brief  density profile from file
 */
template<>
point_t
compute_apm_acc<param::from_file>(body & particle, std::vector<body *> &nbs) {
   using namespace param;
   using namespace density_profiles;

    // Relative density error
    double Pi_a = 0.0, Pi_b = 0.0;
    // Base pressure
    double P0 = apm_base_pressure;

    // Acceleration from artificial pressure
    point_t a_apm = 0.0;

    // this particle (index 'a')
    const double h_a = particle.radius(),
                 rho_a = particle.getDensity();
    const point_t pos_a = particle.coordinates();

    // neighbor particles (index 'b')
    const int n_nb = nbs.size();
    double rho_[n_nb], m_[n_nb], h_[n_nb];
    point_t pos_[n_nb], DiWa_[n_nb];
    
    // compute kernel gradient
    for (int b = 0; b < n_nb; ++b) {
      const body * const nb = nbs[b];
      rho_[b] = nb->getDensity();
      pos_[b] = nb->coordinates();
      h_[b] = nb->radius();
      m_[b] = nb->mass()*(pos_[b] != pos_a);
      const point_t pos_ab = pos_a - pos_[b];
      const double h_ab = 0.5*(h_a + h_[b]);
      DiWa_[b] = kernels::sph_kernel_gradient(pos_ab,h_ab);
    }

    //Invoke Lane-Emden solver to compute density profile
    #ifdef LE_density
    int N_r = lane_emden_radial_N;
    lane_emden::solve(N_r, rad_grid, rho_grid, mass_grid, drhodr_grid,
                      alpha2_grid, dalpha2dr_grid, beta2_grid, dbeta2dr_grid);
    #endif

    //target density profile
    const double r_a = flecsi::magnitude(pos_a);
    double rho_a_target = rho_from_data_grid(r_a);
    
    // compute apm acceleration
    for (int b = 0; b < n_nb; ++b) {
      const double r_b = flecsi::magnitude(pos_[b]);
      double rho_b_target = rho_from_data_grid(r_b);
      Pi_a = std::max(1.0 + (rho_a - rho_a_target)/rho_a_target, 0.1);
      Pi_b = std::max(1.0 + (rho_[b] - rho_b_target)/rho_b_target, 0.1);
      a_apm = -P0*m_[b] * (Pi_a + Pi_b)/(rho_a*rho_[b]) * DiWa_[b];
    }

    return a_apm;

}
#ifdef apm_type
#define sph_compute_apm_acc compute_apm_acc<param::apm_type>
#else
compute_apm_acc_t sph_compute_apm_acc = nullptr;
#endif
/**
 * @brief      APM selector
 */
void
select() {
  using namespace param;
#ifndef apm_type
  switch(apm_type) {
    case(zero_apm):
      sph_compute_apm_acc = compute_apm_acc<zero_apm>;
      break;
    case(kn_ejecta):
      sph_compute_apm_acc = compute_apm_acc<kn_ejecta>;
      break;
    case(sharp_spherical):
      sph_compute_apm_acc = compute_apm_acc<sharp_spherical>;
      break;
    case(from_file):
      sph_compute_apm_acc = compute_apm_acc<from_file>;
      break;
    default:
      log_fatal("Bad APM parameter" << std::endl);
  }
#endif

} // select()

}; // namespace apm
#undef SQ
#undef QU

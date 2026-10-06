/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2017 Triad National Security, LLC
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
 * @file main_driver.cc
 * @author Julien Loiseau
 * @date April 2017
 * @brief Main driver for FleCSPH simulation.
 * The Specialization Driver is normally used to register data and the main
 * code is in the Driver.
 */

#include <iostream>
#include <numeric> // For accumulate

#include <mpi.h>
#include <omp.h>

#undef fmm_order
#include "analysis.h"
#include "bodies_system.h"
#include "default_physics.h"
#include "diagnostic.h"
#include "params.h"
#include "influx.h"

#include "main.h"

#define OUTPUT_ANALYSIS

static std::string output_h5data_file; // = output_h5data_prefix + ".h5part"

void
set_derived_params() {
  using namespace param;

  // set kernel
  kernels::select();

  // set viscosity
  viscosity::select();

  // filenames (this will change for multiple files output)
  std::ostringstream oss;
  oss << output_h5data_prefix << ".h5part";
  output_h5data_file = oss.str();

  // analysis: set output times
  analysis::set_initial_time_iteration();

  // set equation of state
  eos::select();
  
  // set density profile
  density_profiles::select();

  // set external force
  external_force::select(external_force_type);

  if (enable_inflow) {
    // initialize the input flux
    influx::init();
  }

  //TODO : Separate KN relaxation and other apm 
  // set apm select
  apm::select();
   if(do_apm && (apm_type == kn_ejecta)){
      SET_PARAM(sphere_radius, (2.*flow_velocity*kn_ejecta_epoch));
   }
}

///**
// * @brief  invert particles inside the extraction radius from  file
// * @param  ifname - hdf5 file: 3 dimensional array with relevant fields
// *              (density, velocities, temp, etc...)
// *             and spherical grid coordinates
// *         R_ex - const: extraction radius in cm
// */
//void
//spherical_inversion(const char * fileprefix, const double R_ex, const double v_ex, body_system<double, gdimension> & bs) {
//   //step 0: convert body_system into vector of bodies
//   std::vector<body> & particles = bs.getLocalbodies();
//   //step 1: loop through particles
//   for(auto pt : particles){
//      // grabbing particle coords
//      point_t rp = pt.coordinates();
//      //step 2: get relevant field info for particle
//      double vr = density_profiles::Q_ndim_from_data_grid(rp, density_profiles::vr_interp);
//      double vt = density_profiles::Q_ndim_from_data_grid(rp, density_profiles::vt_interp);
//      double vp = density_profiles::Q_ndim_from_data_grid(rp, density_profiles::vp_interp);
//      point_t newvel{vr,vt,vp};
//      pt.setVelocity(newvel);
//      printf()
//      // code for other fields
//      // write interpolated values for other fields
//
//      //step 3: readjust position for correct for real radial velocity
//      // calc spherical coords for particle
//      double x=rp[0], y=rp[1], z=rp[2];
//      double r = sqrt(x*x + y*y + z*z);
//      double theta = atan2(sqrt(x*x+y*y),z);
//      double phi = atan2(y,x);
//      // new radial coord
//      double r_dag = R_ex + ((r-R_ex)/v_ex)*vr;
//
//      //step 4: perform spherical inversion
//      // only radial position changes
//      double r_inv = R_ex*R_ex/r_dag;
//      // write new value to point
//      double x_inv = r_inv*sin(theta)*cos(phi);
//      double y_inv = r_inv*sin(theta)*sin(phi);
//      double z_inv = r_inv*cos(theta);
//      point_t pt_inv{x_inv,y_inv,z_inv};
//      pt.set_coordinates(pt_inv);
//   }
// 
// //step 5: save particles to new output
// 
// // adding tag to file prefix
//    char ofname[128];
//    sprintf(ofname, "%s_inversion", fileprefix);
//    bs.write_bodies(ofname, 0, 0);
// }


int
advance(const std::string& parameter_file) {

  using namespace param;

  int rank;
  int size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // set simulation parameters
  param::mpi_read_params(parameter_file);
  set_derived_params();

  // read input file and initialize equation of state
  body_system<double, gdimension> bs;
  bs.read_bodies(initial_data_prefix, output_h5data_prefix, initial_iteration);
  if (not adaptive_timestep && initial_dt > 0) physics::dt = initial_dt;

  MPI_Barrier(MPI_COMM_WORLD);

  do {
    analysis::screen_output(rank);
    MPI_Barrier(MPI_COMM_WORLD);

    if(physics::iteration == param::initial_iteration) {

      log_one(trace) << "Initial iteration" << std::endl;
      bs.update_iteration();

      // for relaxation phase, reset equation of state to polytropic
      // reset polytropic gamma to 0.99
      if (physics::iteration < relaxation_steps) {
        SET_PARAM(eos_type, eos_polytropic);
        SET_PARAM(poly_gamma, 1.01);
        eos::select();
        body pt0;
        pt0.setDensity(rho_initial);
        pt0.setPressure(pressure_initial);
        eos::compute_entropy(pt0);
        double K = pt0.getEntropy();
        bs.apply_all([&](body & pt) {pt.setEntropy(K);});
        bs.apply_all(eos::compute_pressure);
        bs.apply_all(eos::compute_internal_energy);
        SET_PARAM(relaxation_beta, 
            relaxation_beta*sqrt(pressure_initial/rho_initial)/sphere_radius);
        log_one(info) << "Relaxation beta set to "<< relaxation_beta <<"\n";
      }

      bs.apply_all(eos::compute_entropy);

      if(thermokinetic_formulation) {
        // compute total energy for every particle
        bs.apply_all(physics::set_total_energy);
      }

      if (sph_viscosity != visc_constant) {
        bs.apply_all(viscosity::initialize_alpha);
      }

      log_one(trace) << "compute density pressure cs"<<std::endl;
      bs.apply_in_smoothinglength(physics::compute_density_pressure_soundspeed);
      bs.apply_all(integration::save_velocityhalf);

      if (compute_density_diff_instead_hrate) {
        bs.apply_all(physics::set_density_diff);
      }

      if (sph_viscosity != visc_constant) {
        log_one(trace) << "computing adaptive viscosity" << std::endl;
        bs.apply_in_smoothinglength(viscosity::compute_alpha);
      }

      // compute acceleration
      log_one(trace) << "compute rhs of evolution equations" << std::endl;
      bs.reset_ghosts();
      bs.apply_in_smoothinglength(physics::compute_acceleration);
      if (do_apm) {
        log_one(trace)<<"position is updated via APM criteria"<< std::endl;
        bs.apply_in_smoothinglength(physics::compute_apm_position_correction);
      }
      if (physics::iteration < relaxation_steps) {
        log_one(trace) << "add relaxation terms" << std::endl;
        bs.apply_all(physics::add_drag_acceleration);
        bs.apply_in_smoothinglength(physics::add_short_range_repulsion);
        log_one(trace) << "relaxation terms: done" << std::endl;
      }

      if (adaptive_timestep) {
        // Update timestep in the very beginning
        log_one(trace) << "compute adaptive timestep" << std::endl;
        bs.apply_all(physics::compute_dt);
        bs.get_all(physics::set_adaptive_timestep);
        log_one(trace) << "adaptive timestep: done" << std::endl;
      }

      if (evolve_internal_energy) {
        if (thermokinetic_formulation){
          // compute de/dt
          for (int m=1; m<=pressure_updates_number;++m) { // 1 or 2 passes
            log_one(trace) << "compute dedt: pass " << m  << std::endl;
            bs.apply_in_smoothinglength(physics::compute_dedt);
            if(add_heating_source)
              bs.apply_all(physics::add_heatrate_dedt);
            if (physics::iteration < relaxation_steps)
              bs.apply_all(physics::add_drag_dedt);

            bs.apply_all(physics::recompute_pressure_soundspeed_thermokinetic);
            if (m < pressure_updates_number)
              bs.reset_ghosts(); // skip syncing with the last pass
          }
        }
        else {
          // or compute du/dt
          for (int m=1; m<=pressure_updates_number;++m) { // 1 or 2 passes
            log_one(trace) << "compute dudt: pass " << m  << std::endl;
            bs.apply_in_smoothinglength(physics::compute_dudt);
            if(add_heating_source)
              bs.apply_all(physics::add_heatrate_dudt);
            if (physics::iteration < relaxation_steps)
              bs.apply_all(physics::add_drag_dudt);
            bs.apply_all(physics::recompute_pressure_soundspeed);
            if (m < pressure_updates_number)
              bs.reset_ghosts(); // skip syncing with the last pass
          }
        }
      } // if evolve_internal_energy
      log_one(trace) << "compute initial rhs terms: done" << std::endl;
    }
    else { // not the initial iteration
      log_one(trace) << "leapfrog: kick one" << std::endl;
      if (evolve_internal_energy) {
        if (thermokinetic_formulation)
          bs.apply_all(integration::leapfrog_kick_e);
        else
          bs.apply_all(integration::leapfrog_kick_u);
      }
      bs.apply_all(integration::leapfrog_kick_v);
      bs.apply_all(integration::save_velocityhalf);
      log_one(trace) << "kick one: done" << std::endl;

      log_one(trace) << "leapfrog: drift" << std::endl;
      bs.apply_all(integration::leapfrog_drift);
      log_one(trace) << "drift: done" << std::endl;

      // sync velocities
      bs.update_iteration();
      log_one(trace) << "compute density pressure cs" << std::endl;
      bs.apply_in_smoothinglength(physics::compute_density_pressure_soundspeed);

      if (sph_viscosity != visc_constant) {
        log_one(trace) << "compute adaptive viscosity" << std::endl;
        bs.apply_in_smoothinglength(viscosity::compute_alpha);
      }

      // compute acceleration
      log_one(trace) << "leapfrog: kick two (velocity)" << std::endl;
      bs.reset_ghosts();
      bs.apply_in_smoothinglength(physics::compute_acceleration);
      if (do_apm) {
        log_one(trace)<<"position is updated via APM criteria"<< std::endl;
        bs.apply_in_smoothinglength(physics::compute_apm_position_correction);
      }
      if(physics::iteration < relaxation_steps) {
        bs.apply_all(physics::add_drag_acceleration);
        bs.apply_in_smoothinglength(physics::add_short_range_repulsion);
      }
      bs.apply_all(integration::leapfrog_kick_v);
      log_one(trace) << "kick two (velocity): done" << std::endl;

      // sync velocities: needed for de/dt
      bs.reset_ghosts();

      if (evolve_internal_energy) {
        log_one(trace) << "leapfrog: kick two (energy)" << std::endl;
        if (thermokinetic_formulation) {
          for (int m=1; m<=pressure_updates_number;++m) { // 1 or 2 passes
            log_one(trace) << "compute dedt: pass " << m  << std::endl;
            bs.apply_in_smoothinglength(physics::compute_dedt);
            if(add_heating_source)
              bs.apply_all(physics::add_heatrate_dedt);
            if (physics::iteration < relaxation_steps)
              bs.apply_all(physics::add_drag_dedt);

            bs.apply_all(physics::recompute_pressure_soundspeed_thermokinetic);
            if (m < pressure_updates_number)
              bs.reset_ghosts(); // skip syncing with the last pass
          }

          bs.apply_all(integration::leapfrog_kick_e);
        }
        else {
          for (int m=1; m<=pressure_updates_number;++m) { // 1 or 2 passes
            log_one(trace) << "compute dudt: pass " << m  << std::endl;
            bs.apply_in_smoothinglength(physics::compute_dudt);
            if(add_heating_source)
              bs.apply_all(physics::add_heatrate_dudt);
            if (physics::iteration < relaxation_steps)
              bs.apply_all(physics::add_drag_dudt);

            bs.apply_all(physics::recompute_pressure_soundspeed);
            if (m < pressure_updates_number)
              bs.reset_ghosts(); // skip syncing with the last pass
          }
          bs.apply_all(integration::leapfrog_kick_u);
        }
        log_one(trace) << "kick two (energy): done" << std::endl;
      } // evolve internal energy
    } // not initial iteration

    if(sph_variable_h){
      log_one(trace) << "updating smoothing length" << std::endl;
      bs.get_all(physics::compute_smoothinglength);
      log_one(trace) << ".done" << std::endl;
    }else if(sph_update_uniform_h){
      // The particles moved, compute new smoothing length
      log_one(trace) << "updating smoothing length" << std::endl;
      bs.get_all(physics::compute_average_smoothinglength,bs.getNBodies());
      log_one(trace) << ".done" << std::endl;
    }

    // Compute and output scalar reductions and diagnostic
    // printf("Saving scalar output...\n");
    analysis::scalar_output(bs,rank);
    // printf("Done.\n");
    // printf("Saving hdf5 output...\n");
    analysis::h5data_output(bs, rank);
    // printf("Done.\n");
    // printf("Saving diagnostic output...\n");
    diagnostic::output(bs,rank);
    // printf("Done.\n");

    // Check for nans
    bs.apply_all(physics::check_nans);
    bs.apply_all(physics::check_negativity);

    if(adaptive_timestep) {
      // Update timestep
      log_one(trace) << "compute adaptive timestep" << std::endl;
      bs.apply_all(physics::compute_dt);
      bs.get_all(physics::set_adaptive_timestep);
      log_one(trace) << ".done" << std::endl;
    }
    
    MPI_Barrier(MPI_COMM_WORLD);

    physics::advance_time();

  } while(not physics::termination_criteria());

  return 0; 
} // advance

bool
check_conservation(const std::vector<analysis::e_conservation> & check) {
  return analysis::check_conservation(check);
}

int
main(int argc, char * argv[]) {
  return driver_main(argc, argv);
}

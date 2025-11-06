/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2017 Triad National Security, LLC
 * All rights reserved.
 *~--------------------------------------------------------------------------~*/

/*~--------------------------------------------------------------------------~*
 *
 * /@@@@@@@@ @@           @@@@@@   @@@@@@@@ @@@@@@@  @@      @@
 * /@@///// /@@          @@////@@ @@////// /@@////@@/@@     /@@
 * /@@      /@@  @@@@@  @@    // /@@       /@@   /@@/@@     /@@
 * /@@@@@@@ /@@ @@///@@/@@       /@@@@@@@@@/@@@@@@@ /@@@@@@@@@@
 * /@@////  /@@/@@@@@@@/@@       ////////@@/@@////  /@@//////@@
 * /@@      /@@/@@//// //@@    @@       /@@/@@      /@@     /@@
 * /@@      /@@@//@@@@@@ //@@@@@@  @@@@@@@@ /@@      /@@     /@@
 * //       ///  //////   //////  ////////  //       //      //
 *
 *~--------------------------------------------------------------------------~*/

/**
 * @file analysis.h
 * @author Julien Loiseau
 * @date April 2017
 * @brief Basic analysis
 */

#pragma once

#include "bodies_system.h"
#include "params.h"
#include "wvt.h"
#include <vector>

namespace analysis {

enum e_conservation : size_t {
  MASS = 0,
  ENERGY = 1,
  MOMENTUM = 2,
  ANG_MOMENTUM = 3
};

point_t linear_momentum;
point_t total_ang_mom;
double total_mass;
double total_energy;
double total_kinetic_energy;
double total_internal_energy;
double total_gravitational_energy;
double velocity_part;
double total_density_diff;

/**
 * @brief      Compute the linear momentum
 *
 * @param      bodies  Vector of all the local bodies
 */
void
compute_total_momentum(std::vector<body> & bodies) {
  linear_momentum = {0};
  for(size_t i = 0; i < bodies.size(); ++i) {
    if(bodies[i].type() != NORMAL)
      continue;
    linear_momentum += bodies[i].mass() * bodies[i].getVelocity();
  }
  mpi_utils::reduce_sum(linear_momentum);
}

/**
 * @brief      Compute the total mass
 *
 * @param      bodies  Vector of all the local bodies
 */
void
compute_total_mass(std::vector<body> & bodies) {
  total_mass = 0.;
  for(size_t i = 0; i < bodies.size(); ++i) {
    if(bodies[i].type() != NORMAL)
      continue;
    total_mass += bodies[i].mass();
  }
  mpi_utils::reduce_sum(total_mass);
}

/**
 * @brief      Compute the total energy (internal + kinetic)
 *
 * @param      bodies  Vector of all the local bodies
 */
void
compute_total_energy(std::vector<body> & bodies) {
  using namespace param;

  total_energy = 0.;
  if (thermokinetic_formulation) {
    for(size_t i = 0 ; i < bodies.size(); ++i) {
      body & pt = bodies[i];
      if(pt.type() != NORMAL)  continue;
      total_energy += bodies[i].mass()*bodies[i].getTotalenergy();
    }
  }
  else {
    for(size_t i = 0 ; i < bodies.size(); ++i){
      body & pt = bodies[i];
      if(pt.type() != NORMAL)  continue;
      const point_t
          pos = pt.coordinates(),
          vel = pt.getVelocity();
      const double
          m = pt.mass(),
          eint = pt.getInternalenergy(),
          epot = external_force::potential(pos),
          ekin = .5*flecsi::dot(vel,vel);
      total_energy += m*(ekin + eint + epot);
    }
  }
  if(enable_fmm) {
    for(size_t i = 0; i < bodies.size(); ++i) {
      body & pt = bodies[i];
      if(pt.type() != NORMAL)  continue;
      // factor of 0.5 takes care of double-counting
      total_energy += 0.5*pt.getGPotential()*pt.mass();
    }
  }
  mpi_utils::reduce_sum(total_energy);
}

/**
 * @brief      Sum up total kinetic energy
 *
 * @param      bodies  Vector of all the local bodies
 */
void
compute_total_kinetic_energy(std::vector<body>& bodies) {
  using namespace param;

  total_kinetic_energy = 0.;
  for(size_t i = 0 ; i < bodies.size(); ++i){
    body & pt = bodies[i];
    if(pt.type() != NORMAL)  continue;
    const double m = pt.mass();
    const point_t vel = pt.getVelocity();
    total_kinetic_energy += .5*m*flecsi::dot(vel,vel);
  }
  mpi_utils::reduce_sum(total_kinetic_energy);
}

/**
 * @brief      Sum up internal energy
 * @param      bodies  Vector of all the local bodies
 */
void
compute_total_internal_energy(std::vector<body> & bodies) {
  using namespace param;

  total_internal_energy = 0.;
  for(size_t i = 0; i < bodies.size(); ++i) {
    if(bodies[i].type() != NORMAL)
      continue;
    total_internal_energy += bodies[i].mass() * bodies[i].getInternalenergy();
  }
  mpi_utils::reduce_sum(total_internal_energy);
}

/**
 * @brief      Sum up gravitational energy
 * @param      bodies  Vector of all the local bodies
 */
void
compute_total_gravitational_energy(std::vector<body> & bodies) {
  using namespace param;

  total_gravitational_energy = 0.;
  for(size_t i = 0; i < bodies.size(); ++i) {
    if(bodies[i].type() != NORMAL)
      continue;
    total_gravitational_energy += bodies[i].getGPotential()*bodies[i].mass();
  }
  // account for double-counting
  total_gravitational_energy *= 0.5;
  mpi_utils::reduce_sum(total_gravitational_energy);
}

/**
 * @brief      Compute total angular momentum
 *
 * @param      bodies  Vector of all the local bodies
 */
void
compute_total_ang_mom(std::vector<body> & bodies) {
  total_ang_mom = {0};
  if constexpr(gdimension == 2) {
    for(size_t i = 0; i < bodies.size(); ++i) {
      if(bodies[i].type() != NORMAL)
        continue;
      const double m = bodies[i].mass();
      const point_t v = bodies[i].getVelocity();
      const point_t r = bodies[i].coordinates();
      total_ang_mom[0] += m*(r[0]*v[1] - r[1]*v[0]);
    }
    mpi_utils::reduce_sum(total_ang_mom);
  }
  if constexpr(gdimension == 3) {
    for(size_t i = 0; i < bodies.size(); ++i) {
      if(bodies[i].type() != NORMAL)
        continue;
      const double m = bodies[i].mass();
      const point_t v = bodies[i].getVelocity();
      const point_t r = bodies[i].coordinates();
      total_ang_mom[0] += m*(r[1]*v[2] - r[2]*v[1]);
      total_ang_mom[1] += m*(r[2]*v[0] - r[0]*v[2]);
      total_ang_mom[2] += m*(r[0]*v[1] - r[1]*v[0]);
    }
    mpi_utils::reduce_sum(total_ang_mom);
  }
}

/**
 * @brief      Compute the sum of density difference between particles and given "unfolded" density
 *
 * @param      bodies  Vector of all the local bodies
 */
void
compute_cumulative_density_diff(std::vector<body> & bodies) {
  total_density_diff = 0.;
  int out_count = 0;
  for(size_t i = 0; i < bodies.size(); ++i) {
    if(bodies[i].type() != NORMAL)
      continue;
    const point_t pos = bodies[i].coordinates();
    const double rho_prof = density_profiles::density_ndim(pos);
    const double avg_rho_err = abs(bodies[i].getDensity() - rho_prof)/rho_prof;
    if (sqrt(pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]) < param::sphere_radius) { 
      if (avg_rho_err > 1e4) {
        out_count++;
        printf("Found outlier particle #%d. Omitting from cumulative err... \n", out_count);
      }
      else {
        total_density_diff += avg_rho_err/(bodies.size()-out_count);
      }
    }
  }

  mpi_utils::reduce_sum(total_density_diff);
}


/**
 * @brief Initialize output times
 */
void
set_initial_time_iteration() {
  using namespace param;
  using namespace physics;

  // iteration and time
  iteration = initial_iteration;
  totaltime = initial_time;
  totaltime_prev = initial_time;
  dt = initial_dt;
  dt_saved = 0.0;

  if (out_screen_dt > 0.0) { // set next screen output time
    t_screen_output = out_screen_dt*((int64_t)(totaltime/out_screen_dt));
    if (t_screen_output < totaltime)
      t_screen_output = totaltime;
  }

  if (out_scalar_dt > 0.0) { // set next scalar output time
    t_scalar_output = out_scalar_dt*((int64_t)(totaltime/out_scalar_dt));
    if (t_scalar_output < totaltime)
      t_scalar_output = totaltime;
  }

  if (out_h5data_dt > 0.0) { // set next h5data output time
    t_h5data_output = out_h5data_dt*((int64_t)(totaltime/out_h5data_dt));
    if (t_h5data_output < totaltime)
      t_h5data_output = totaltime;
  }

}


/**
 * @brief Rolling screen output
 */
void
screen_output(int rank) {
  using namespace param;
  static int count = 0;
  const int screen_length = 40;

  if (out_screen_dt > 0.0) { // output by time
    if (physics::totaltime < physics::t_screen_output) {
      return;
    }
  }
  else { // output by iteration
    if (out_screen_every <= 0)
      return;
    if (physics::iteration % out_screen_every != 0)
      return;
  }
  (++count - 1) % screen_length ||
    log_one(info) << "#-- iteration:               time:" << std::endl;
  log_one(info) << std::setw(14) << physics::iteration << std::setw(20)
                << std::scientific << std::setprecision(12)
                << physics::totaltime << std::endl;
}

/**
 * @brief Periodic file output
 *
 * Outputs all scalar reductions as single formatted line to a file. When
 * creating the file, writes a header to indicate which quantities are
 * output in which column (because their order and quantity may change
 * between revisions)
 * E.g.:
 * -- >> example output file >> -----------------------------------------
 * # Scalar reductions:
 * # 1:iteration 2:time 3:energy 4:mom_x 5:mom_y 6:mom_z
 * 0  0.0   1.0000000e+03  0.00000000e+00  0.00000000e+00  0.00000000e+00
 * 10 0.1   1.0000000e+03  0.00000000e+00  0.00000000e+00  0.00000000e+00
 * 20 0.2   1.0000000e+03  0.00000000e+00  0.00000000e+00  0.00000000e+00
 * 30 0.3   1.0000000e+03  0.00000000e+00  0.00000000e+00  0.00000000e+00
 * ...
 * -- << end output file <<<< -------------------------------------------
 */
void
scalar_output(body_system<double, gdimension> & bs, const int rank) {
  static bool first_time = true;
  if (param::out_scalar_dt > 0.0) { // output by time
    if (physics::totaltime < physics::t_scalar_output) {
      return;
    }
  }
  else { // output by iteration
    if (param::out_scalar_every <= 0)
      return;
    if (physics::iteration % param::out_scalar_every != 0)
      return;
  }

  // // TODO: maybe put this into termination_criteria
  // if(wvt_basic::wvt_converged == false &&
  //    physics::iteration != (param::final_iteration + param::wvt_cool_down))
  //   return;

  // compute reductions
  bs.get_all(compute_total_momentum);
  bs.get_all(compute_total_mass);
  bs.get_all(compute_total_energy);
  bs.get_all(compute_total_kinetic_energy);
  bs.get_all(compute_total_internal_energy);
  bs.get_all(compute_total_gravitational_energy);
  bs.get_all(compute_total_ang_mom);
  bs.get_all(compute_cumulative_density_diff);

  // output only from rank #0
  if(rank != 0)
    return;
  const char * filename = "scalar_reductions.dat";

  if(first_time) {
    // Generate and output the header
    std::ostringstream oss_header;
    switch(gdimension) {
      case 1:
        oss_header << "# Scalar reductions: " << std::endl
                   << "# 1:iteration 2:time 3:timestep 4:total_mass "
                   << " 5:total_energy 6:kinetic_energy 7:internal_energy"
                   << " 8:mom_x" << std::endl;
        break;

      case 2:
        oss_header
          << "# Scalar reductions: " << std::endl
          << "# 1:iteration 2:time 3:timestep 4:total_mass 5:total_energy"
          << " 6:kinetic_energy 7:internal_energy" << std::endl
          << "# 8:mom_x 9:mom_y 10:and_mom_z" << std::endl;
        break;

      case 3:
      default:
        if (param::enable_fmm and not(param::evolve_internal_energy)) {
          if (param::compute_density_diff_instead_hrate)
            oss_header
              << "# Scalar reductions: " << std::endl
              << "# 1:iteration 2:time 3:timestep 4:total_mass 5:total_energy"
              << " 6:kinetic_energy 7:gravitational_energy " << std::endl
              << "# 8:fractional_density_err:" << std::endl
              << "# 9:mom_x 10:mom_y 11:mom_z "
              << "12:ang_mom_x 13:ang_mom_y 14:ang_mom_z" << std::endl
              << "# 15: com_x 16: com_y 17: com_z" << std::endl;
           else
             oss_header
               << "# Scalar reductions: " << std::endl
               << "# 1:iteration 2:time 3:timestep 4:total_mass 5:total_energy"
               << " 6:kinetic_energy 7:gravitational_energy " << std::endl
               << "# 8:mom_x 9:mom_y 10:mom_z "
               << "11:ang_mom_x 12:ang_mom_y 13:ang_mom_z" << std::endl
               << "# 14: com_x 15: com_y 16: com_z" << std::endl;
         }

        else {
          if (param::compute_density_diff_instead_hrate)
            oss_header
              << "# Scalar reductions: " << std::endl
              << "# 1:iteration 2:time 3:timestep 4:total_mass 5:total_energy"
              << " 6:kinetic_energy 7:gravitational_energy " << std::endl
              << "# 8:fractional_density_err:" << std::endl
              << "# 9:mom_x 10:mom_y 11:mom_z "
              << "12:ang_mom_x 13:ang_mom_y 14:ang_mom_z" << std::endl
              << "# 15: com_x 16: com_y 17: com_z" << std::endl; 
          else
            oss_header
              << "# Scalar reductions: " << std::endl
              << "# 1:iteration 2:time 3:timestep 4:total_mass 5:total_energy"
              << " 6:kinetic_energy 7:internal_energy " << std::endl
              << "# 8:mom_x 9:mom_y 10:mom_z "
              << "11:ang_mom_x 12:ang_mom_y 13:ang_mom_z" << std::endl
              << "# 14: com_x 15: com_y 16: com_z" << std::endl;
        }

    }

    std::ofstream out(filename);
    out << oss_header.str();
    out.close();
    first_time = false;
  }

  std::ostringstream oss_data;
  oss_data << std::setw(14) << physics::iteration << std::setw(20)
           << std::scientific << std::setprecision(12) << physics::totaltime
           << std::setw(20) << physics::dt << " " << total_mass << " "
           << total_energy << " " << total_kinetic_energy << " ";
  // if internal energy is not evolved, output gravitational energy
  if (param::enable_fmm and not(param::evolve_internal_energy))
      oss_data << total_gravitational_energy << " ";
  else
      oss_data << total_internal_energy << " ";

  if (param::compute_density_diff_instead_hrate)
      oss_data << total_density_diff << " ";

  for(unsigned short int k = 0; k < gdimension; ++k)
    oss_data << " " << linear_momentum[k];

  if(gdimension == 2)
    oss_data << " " << total_ang_mom[0];

  if(gdimension == 3) {
    for(unsigned short k = 0; k < gdimension; ++k)
      oss_data << " " << total_ang_mom[k];
    for(unsigned short k = 0; k < gdimension; ++k)
      oss_data << " " << bs.tree()->root_node()->coordinates()[k];
  }

  oss_data << std::endl;

  // Open file in append mode
  std::ofstream out(filename, std::ios_base::app);
  out << oss_data.str();
  out.close();

} // scalar output

/**
 * @brief Periodic output to an h5part file
 */
void
h5data_output(body_system<double, gdimension> & bs, const int rank) {
  using namespace param;
  using namespace physics;

  if (out_h5data_dt > 0.0) { // output by time
    if (totaltime < t_h5data_output) {
      return;
    }
  }
  else { // output by iteration
    if (out_h5data_every <= 0)
      return;
    if (iteration % out_h5data_every != 0)
      return;
  }
  bs.write_bodies(output_h5data_prefix, iteration, totaltime);
} // h5data_output

bool
check_conservation(const std::vector<e_conservation> & check) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  log_one(info) << "Checking conservation of: ";
  for(auto c : check) {
    switch(c) {
      case MASS:
        log_one(info) << " MASS ";
        break;
      case ENERGY:
        log_one(info) << " ENERGY ";
        break;
      case MOMENTUM:
        log_one(info) << " MOMENTUM ";
        break;
      case ANG_MOMENTUM:
        log_one(info) << " ANG_MOMENTUM ";
        break;
      default:
        break;
    } // switch
  }
  // Open file and check
  double tol = 1.0e-16;
  size_t iteration;
  double time, timestep, mass, energy, momentum[gdimension], ang_mom[3];
  double base_energy, base_mass, base_ang_mom[3], base_momentum[gdimension];
  std::ifstream inFile;
  inFile.open("scalar_reductions.dat");
  if(!inFile) {
    std::cerr << "Unable to open file scalar_reductions.dat";
    exit(1); // call system to stop
  }

  // Read the first line
  inFile >> iteration >> time >> timestep >> base_mass >> base_energy;
  for(size_t i = 0; i < gdimension; ++i) {
    inFile >> base_momentum[i];
  }
  if(gdimension == 2) {
    inFile >> base_ang_mom[2];
  }
  if(gdimension == 3) {
    inFile >> base_ang_mom[0] >> base_ang_mom[1] >> base_ang_mom[2];
  }

  while(inFile.good()) {
    inFile >> iteration >> time >> timestep >> mass >> energy;
    for(size_t i = 0; i < gdimension; ++i) {
      inFile >> momentum[i];
    }
    if(gdimension == 2) {
      inFile >> ang_mom[2];
    }
    if(gdimension == 3) {
      inFile >> ang_mom[0] >> ang_mom[1] >> ang_mom[2];
    }
    // Check the quantities
    // Check only the required ones
    for(auto c : check) {
      switch(c) {
        case MASS:
          if(!(abs(base_mass - mass) < tol)) {
            std::cerr << "Mass is not conserved" << std::endl;
            return false;
          }
          break;
        case ENERGY:
          if(!(abs(base_energy - energy) < tol)) {
            std::cerr << "Energy is not conserved" << std::endl;
            return false;
          }
          break;
        case MOMENTUM:
          if(!(abs(base_momentum - momentum) < tol)) {
            std::cerr << "Momentum is not conserved" << std::endl;
            return false;
          }
          break;
        case ANG_MOMENTUM:
          if(!(abs(base_ang_mom - ang_mom) < tol)) {
            std::cerr << "Angular Momentum is not conserved" << std::endl;
            return false;
          }
          break;
        default:
          break;
      } // switch
    } // for
  } // while
  // Close the file
  inFile.close();
  return true;
} // conservation check

}; // namespace analysis

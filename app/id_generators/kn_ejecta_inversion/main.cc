/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2025 Triad National Security, LLC
 * All rights reserved.
 *~--------------------------------------------------------------------------~*/

#include <algorithm>
#include <cassert>
#include <iostream>
#include <math.h>
#include <random>

#include "density_profiles.h"
#include "io.h"
#include "kernels.h"
#include "lattice.h"
#include "params.h"
#include "sedov.h"
#include "user.h"
#include "eos.h"
using namespace io;
#include "bodies_system.h"
#include "influx.h"
#include "analysis.h"

#define SQ(x) ((x) * (x))
#define CU(x) ((x) * (x) * (x))
#define QU(x) ((x) * (x) * (x) * (x))

/*
  Create a 3D compactified distribution from an existing particle
  input data.
*/

//
// help message
//
void
print_usage() {
  std::cout << "\"Invert\" particle configuration to create a compactified\n"
            << "initial data (such as the one made by Brendan)\n"
            << "Usage: ./kn_ejecta_3d_generator <parameter-file.par>"
            << std::endl;
}

//
// derived parameters
//
static double timestep = 1.0; // recommended timestep
static double total_mass = 1.; // total mass of the fluid
static double mass_particle = 1.; // mass of an individual particle
static point_t bbox_max, bbox_min; // bounding box of the domain
static char initial_data_file[256]; // = initial_data_prefix[_XXXXX].h5part"
static const double domain_type_sphere = 1;
static double rho_c = 1.; // central density

void
set_derived_params() {
  using namespace param;

  if constexpr (gdimension != 3) {
    log_one(error) << "this initial data only works in 3D" << std::endl;
    MPI_Abort(MPI_COMM_WORLD, -1);
  }

  // set kernel
  kernels::select();

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

}

int
main(int argc, char * argv[]) {
  using namespace param;
  using namespace influx;

  // check options list: exactly one option is allowed
  if(argc != 2) {
    std::cerr << "ERROR: parameter file not specified!" << std::endl;
    print_usage();
    exit(0);
  }

  // launch MPI
  int rank, size;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  assert(size == 1); // parallel ID generator not implemented yet
  log_set_output_rank(0);

  // set simulation parameters
  param::mpi_read_params(argv[1]);
  set_derived_params();
  body_system<double, gdimension> bs;
  assert(modify_initial_data); // this is a modifier, so "modify_initial_data" must be "yes"
  bs.read_bodies(initial_data_prefix, "", initial_iteration);
  auto & particles = bs.getLocalbodies();
  std::vector<body> new_particles{};


  //printf("initial iteration passed to id generator is: %5d\n",initial_iteration);
  // go through the particles and set interpolated quantities
  for(auto & pt : particles) {
     // grabbing particle coords
     point_t rp = pt.coordinates();
     //step 2: get relevant field info for particle
     double vr = density_profiles::Q_ndim_from_data_grid(rp, density_profiles::vr_interp)*C_LIGHT_CGS;
     if (std::abs(vr)<1e-15*C_LIGHT_CGS) continue;

     double vt = density_profiles::Q_ndim_from_data_grid(rp, density_profiles::vt_interp)*C_LIGHT_CGS;
     double vp = density_profiles::Q_ndim_from_data_grid(rp, density_profiles::vp_interp)*C_LIGHT_CGS;
     point_t newvel{vr,vt,vp};
     pt.setVelocity(newvel);
     // printf("the expansion radial velocity is: %12.5e \n", flow_velocity);
     // printf("the correct radial velocity is: %12.5e \n", vr);

     // code for other fields
     // write interpolated values for other fields

     //step 3: readjust position for correct for real radial velocity
     // calc spherical coords for particle
     double x=rp[0], y=rp[1], z=rp[2];
     double r = sqrt(x*x + y*y + z*z);
    //  printf("the radius before inversion: %12.5e \n",r);
     double theta = atan2(sqrt(x*x+y*y),z);
     double phi = atan2(y,x);
     // new radial coord
     double r_dag = sphere_radius + ((r-sphere_radius)/(flow_velocity*C_LIGHT_CGS))*vr;
     // printf("the adjusted radius before inversion: %12.5e \n",r_dag);

     //step 4: perform spherical inversion
     // only radial position changes
     double r_inv = sphere_radius*sphere_radius/r_dag;
     if (r_inv > sphere_radius) continue;
     // printf("the radius after inversion: %12.5e \n",r_inv);
     // write new value to point
     double x_inv = r_inv*sin(theta)*cos(phi);
     double y_inv = r_inv*sin(theta)*sin(phi);
     double z_inv = r_inv*cos(theta);
     point_t pt_inv{x_inv,y_inv,z_inv};
     // printf("the inverted point is: %12.5e %12.5e %12.5e \n", pt_inv[0], pt_inv[1], pt_inv[2]);
     pt.set_coordinates(pt_inv);
     point_t rp_inv = pt.coordinates();
     // printf("the new coordinates stored in the particle are: %12.5e %12.5e %12.5e \n", rp_inv[0], rp_inv[1], rp_inv[2]);

     if (eos_type == eos_polytropic) {
        pt.setDensity(rho_initial);
        pt.setPressure(pressure_initial);
        eos::compute_entropy(pt);
     }

     pt.set_state(INACTIVE);
     pt.setDensity(density_profiles::Q_ndim_from_data_grid(rp, density_profiles::rho_interp));

     // making sure the internal energy is specific internal energy (per mass):
     pt.setInternalenergy(density_profiles::Q_ndim_from_data_grid(rp, density_profiles::eps_interp));
     pt.setAbar(initial_abar);
     pt.setElectronfraction(initial_zbar/initial_abar);
     eos::compute_pressure(pt);
     eos::compute_temperature(pt);

     if (eos_type == eos_polytropic) eos::compute_internal_energy(pt);
     
     //// pt.setElectronfraction(gp.ye);
     // pt.setPressure(gp.pres);
     // pt.setTemperature(gp.temp);

     // eos::compute_entropy(pt);
     new_particles.push_back(pt);
  }

  SET_PARAM(nparticles, new_particles.size());
  log_one(info) << "Number of particles: " << nparticles << std::endl;

  // write the file; iteration for initial data MUST BE zero!!
  //bs.write_bodies(output_h5data_prefix, 0, 0.0);
  // remove the previous file

  char output_h5data_file[256];
  physics::totaltime = param::initial_time;
  sprintf(output_h5data_file, "%s.h5part", output_h5data_prefix);
  remove(output_h5data_file);
  io::outputDataHDF5(new_particles, output_h5data_prefix, 0, 0.0);
  MPI_Finalize();
  return 0;
}

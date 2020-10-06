/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2017 Triad National Security, LLC
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

#define SQ(x) ((x) * (x))
#define CU(x) ((x) * (x) * (x))

/*
Sets up spherically-symmetric analytic density profile for kilonova ejecta:

  rho (r, t) = rho_0 (t/t0)^-3 (1 - r^2/(v_max*t)^2)^3

Ejecta profile is fully determined by the following parameters:

 - kn_ejecta_mass: total mass of the ejecta (in Msun)
 - kn_ejecta_vmed: median velocity (in clight)
 - kn_ejecta_temp: temperature at the centre

Reference: Wollaeger et al. (2018), arXiv:1705.07084, Section 2.2.1
*/

//
// help message
//
void
print_usage() {
  std::cout << "Initial data generator for analytic kilonova ejecta\n" 
            << "Usage: ./kn_ejecta_generator <parameter-file.par>" 
            << std::endl;
}

//
// derived parameters
//
static double timestep = 1.0; // Recommended timestep
static double total_mass = 1.; // total mass of the fluid
static double mass_particle = 1.; // mass of an individual particle
static point_t bbox_max, bbox_min; // bounding box of the domain
static char initial_data_file[256]; // = initial_data_prefix[_XXXXX].h5part"
static const double domain_type_sphere = 1;

void
set_derived_params() {
  using namespace param;

  if constexpr (gdimension != 3) {
    log_one(error) << "this initial data only works in 3D" << std::endl;
    MPI_Abort(MPI_COMM_WORLD, -1);
  }

  eos::select();
  density_profiles::select();
  particle_lattice::select();

  // The value for constant timestep
  timestep = initial_dt;

  // Bounding box of the domain
  bbox_min = -sphere_radius;
  bbox_max = sphere_radius;

  // particle separation
  SET_PARAM(sph_separation, (2. * sphere_radius / (lattice_nx - 1)));

  // Count number of particles
  int64_t tparticles = particle_lattice::count(
    lattice_type, domain_type_sphere, bbox_min, bbox_max, sph_separation, 0);
  SET_PARAM(nparticles, tparticles);

  // total mass: normalize mass such that central density is rho_initial
  // TODO: instead, derive rho_initial from kn_ejecta_mass
  total_mass = rho_initial * CU(sphere_radius) /
               density_profiles::spherical_density_profile(0.0);

  // single particle mass
  assert(equal_mass);
  mass_particle = total_mass / nparticles;

  // set kernel
  kernels::select();

  // smoothing length
  const double sph_h = sph_eta * kernels::kernel_width *
                       pow(mass_particle / rho_initial, 1. / gdimension);
  SET_PARAM(sph_smoothing_length, sph_h);

  // Filename to be generated
  bool input_single_file = H5P_fileExists(initial_data_prefix);
  if(input_single_file or initial_iteration == 0)
    sprintf(initial_data_file, "%s.h5part", initial_data_prefix);
  else {

    // find the file with initial_iteration
    int step =
      H5P_findIterationSnapshot(initial_data_prefix, param::initial_iteration);
    // file doesn't exist: complain and exit
    if(step < 0) {
      log_one(error) << "Cannot find iteration " << param::initial_iteration
                     << " in prefix " << initial_data_prefix << std::endl;
      exit(MPI_Barrier(MPI_COMM_WORLD) && MPI_Finalize());
    }
    sprintf(initial_data_file, "%s_%05d.h5part", initial_data_prefix, step);
  }
}

int
main(int argc, char * argv[]) {
  using namespace param;

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
  if(modify_initial_data) {
    bs.read_bodies(initial_data_prefix, "", initial_iteration);
    SET_PARAM(nparticles, bs.getNBodies());
  }
  else {
    bs.getLocalbodies().clear();
    bs.getLocalbodies().resize(nparticles);
  }
  auto & bodies = bs.getLocalbodies();

  // Declare coordinate arrays
  double * x = new double[nparticles]();
  double * y = new double[nparticles]();
  double * z = new double[nparticles]();

  if(not modify_initial_data) {
    // Generate the lattice
    auto _np = particle_lattice::generate(lattice_type, domain_type_sphere, 
        bbox_min, bbox_max, sph_separation, 0, x, y, z);
    assert(nparticles == _np);

    for(int64_t a = 0L; a < nparticles; ++a) {
      body & particle = bodies[a];
      point_t pos = {x[a], y[a], z[a]};
      particle.set_coordinates(pos);
    }
  }

  // Assign density, pressure and specific internal energy to particles,
  // including the particles in the blast zone
  const double rho0 = density_profiles::spherical_density_profile(0);

  // For given initial pressure and density, compute adiabatic invariant;
  // this adiabatic invariant is used in the loop below to set up all
  // other thermodynamic quantities ("constant entropy" setup).
  body pt0;
  pt0.setPressure(pressure_initial);
  pt0.setDensity(rho_initial);
  pt0.setAbar(initial_abar);
  pt0.setElectronfraction(initial_zbar/initial_abar);
  pt0.setTemperature(initial_temp);
  eos::compute_entropy(pt0);
  double K0 = pt0.getEntropy();

  // Main loop: assign quantities on particles
  std::default_random_engine generator;
  for(int64_t a = 0; a < nparticles; ++a) {
    body & particle = bodies[a];

    // zero velocity for this test
    point_t zero = 0;
    particle.setVelocity(zero);
    particle.setAcceleration(zero);

    // radial distance from the origin
    point_t rp(particle.coordinates());
    double r = magnitude(rp);

    // set density, particle mass, smoothing length and id
    double rho_a, m_a, h_a;
    if(modify_initial_data) {
      rho_a = particle.getDensity();
      m_a = particle.mass();
      h_a = sph_eta*kernels::kernel_width*pow(m_a/rho_a, 1./gdimension);
    }
    else {
      rho_a = rho_initial / rho0 // renormalize density profile
              * density_profiles::spherical_density_profile(r / sphere_radius);
      m_a = mass_particle;
      h_a = sph_eta * kernels::kernel_width *
            pow(mass_particle / rho_a, 1. / gdimension);
      particle.setDensity(rho_a);
      particle.set_mass(m_a);
      particle.set_radius(h_a);
      particle.set_id(a);
    }

    if(lattice_perturbation_amplitude > 0.0) {
      // add lattice perturbation
      std::normal_distribution<double> distribution(
        0., h_a * lattice_perturbation_amplitude);
      for(unsigned short k = 0; k < gdimension; ++k) {
        rp[k] += distribution(generator);
      }
      particle.set_coordinates(rp);
    }

    // set uniform composition
    particle.setAbar(initial_abar);
    particle.setElectronfraction(initial_zbar/initial_abar);
    particle.setTemperature(initial_temp);

    // set internal energy
    particle.setEntropy(K0);
    eos::compute_internal_energy(particle);

    // set pressure (a function of density and internal energy)
    eos::compute_pressure(particle);

    // set timestep
    particle.setDt(initial_dt);
  }

  log_one(info) << "Number of particles: " << nparticles << "\n"
                << "Total mass:          " << total_mass << std::endl;

  // remove the previous file
  remove(initial_data_file);
  delete[] x;
  delete[] y;
  delete[] z;

  // write the file; iteration for initial data MUST BE zero!!
  bs.write_bodies(initial_data_prefix, 0, 0.0);
  MPI_Finalize();
  return 0;
}

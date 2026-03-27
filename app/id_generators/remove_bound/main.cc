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

/*
  Removes bound particles from the input file and writes the result to output.
*/

//
// help message
//
void
print_usage() {
  std::cout << "Remove bound particles\n"
            << "Usage: ./remove_bound_3d_generator <parameter-file.par>\n";
}

//
// derived parameters
//
static double timestep = 1.0;       // Recommended timestep
static double total_mass = 1.;      // total mass of the fluid
static double mass_particle = 1.;   // mass of an individual particle
static point_t bbox_max, bbox_min;  // bounding box of the domain
static char initial_data_file[256]; // = initial_data_prefix[_XXXXX].h5part"

void
set_derived_params() {
  using namespace param;

  // set external force
  external_force::select(external_force_type);

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
  bs.read_bodies(initial_data_prefix, "", initial_iteration);
  SET_PARAM(nparticles, bs.getNBodies());
  auto & bodies = bs.getLocalbodies();

  // main loop
  int64_t n_removed = 0L;
  for(auto pt = bodies.begin(); pt < bodies.end();) {
    const point_t vel = pt->getVelocity(), pos = pt->coordinates();
    double etot = .5*flecsph::dot(vel,vel) + external_force::potential(pos);
    if (etot < 0) {
      pt = bodies.erase(pt);
      ++n_removed;
    }
    else {
      ++pt;
    }
  }
  log_one(info) << "Removed " << n_removed << " bound particles.\n";
  log_one(info) << "New number of particles: " << bodies.size() << "\n";

  // output
  remove(output_h5data_prefix);
  io::outputDataHDF5(bodies, output_h5data_prefix, 0, 0.);

  // cleanup
  MPI_Finalize();
  return 0;
}

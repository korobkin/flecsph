/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2017 Triad National Security, LLC
 * All rights reserved.
 *~--------------------------------------------------------------------------~*/

#include <algorithm>
#include <cassert>
#include <iostream>
#include <math.h>

#include "io.h"
#include "kernels.h"
#include "lattice.h"
#include "params.h"
#include "KH.h"
#include "user.h"

using namespace io;
#include "bodies_system.h"


//
// help message
//
void
print_usage() {
  log_one(warn) << "Initial data generator for KH test in" << gdimension << "D"
                << std::endl
                << "Usage: ./KD_XD_generator <parameter-file.par>" << std::endl;
}

//
// derived parameters
//
static double rho_m, rho_t; // densities
static double vx_m, vx_t; // velocities
static double pressure_m, pressure_t; // pressures
static char initial_data_file[256]; // = initial_data_prefix[_XXXXX].h5part"

// geometric extents of the three regions: top, middle and bottom
static point_t tbox_min, tbox_max;
static point_t mbox_min, mbox_max;
static point_t bbox_min, bbox_max;

// lattice spacing
// capital letters are for the full-periodicity lattice thickness
static double dx_m, dy_m, dY_m, dz_m, dZ_m; // medium layer
static double dx_t, dy_t, dY_t, dz_t, dZ_t; // medium layer

// width and height of the middle layer, top layer and gap between them
static double w_m, w_t, gap;
static double h_m, h_t;

static int64_t np_middle = 0; // number of particles in the middle block
static int64_t np_top = 0; // number of particles in the top block
static int64_t np_bottom = 0; // number of particles in the bottom block
static double sph_sep_t = 0; // particle separation in top or bottom blocks
static double pmass = 0; // particle mass in the middle block
static double pmass_t = 0; // particle mass in top or bottom blocks

void
set_derived_params() {
  using namespace std;
  using namespace param;

  // support for only equal-mass configurations for now
  if(not equal_mass) {
    log_one(error) << "Only equal-mass configurations are implemented"
                   << std::endl;
    MPI_Finalize();
    exit(0);
  }

  // domain must be rectangular
  assert(domain_type == 0);

  mbox_max[0] = bbox_max[0] = tbox_max[0] = box_length / 2.;
  mbox_min[0] = bbox_min[0] = tbox_min[0] = -box_length / 2.;

  mbox_max[1] = box_width / 4.;
  mbox_min[1] = -box_width / 4.;
  bbox_min[1] = -box_width / 2.;
  bbox_max[1] = -box_width / 4.;
  tbox_min[1] = box_width / 4.;
  tbox_max[1] = box_width / 2.;

  if(gdimension == 3) {
    mbox_max[2] = bbox_max[2] = tbox_max[2] = box_height / 2.;
    mbox_min[2] = bbox_min[2] = tbox_min[2] = -box_height / 2.;
  }

  // set physical parameters
  // --  in the top and bottom boxes (tbox and bbox):
  rho_t = rho_initial;
  pressure_t = pressure_initial;
  vx_t = -flow_velocity / 2.0;
  // -- in the middle box
  rho_m = rho_t * density_ratio; // 2.0 by default
  pressure_m = pressure_t; // pressures must be equal in KH test
  vx_m = flow_velocity / 2.0;

  // file to be generated
  bool input_single_file = H5P_fileExists(initial_data_prefix);
  if(input_single_file or initial_iteration == 0) {
    sprintf(initial_data_file, "%s.h5part", initial_data_prefix);
  }
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

  // select particle lattice and kernel function
  particle_lattice::select();
  kernels::select();

  // particle mass and spacing
  SET_PARAM(sph_separation, box_length / lattice_nx);
  if(gdimension == 3) {
    pmass = rho_m * sph_separation * sph_separation * sph_separation;
    if(lattice_type == 1 or lattice_type == 2)
      pmass *= 1. / sqrt(2.);
    sph_sep_t = sph_separation * cbrt(density_ratio);
  }
  if(gdimension == 2) {
    pmass = rho_m * sph_separation * sph_separation;
    if(lattice_type == 1 or lattice_type == 2)
      pmass *= sqrt(3) / 2;
    sph_sep_t = sph_separation * sqrt(density_ratio);
  }
  pmass_t = pmass;

  // lattice spacing
  dx_m = dy_m = dY_m = dz_m = dZ_m = sph_separation;
  dx_t = dy_t = dY_t = dz_t = dZ_t = sph_sep_t;
  if(lattice_type == 0) {
    log_one(info) << "Lattice: rectangular, resolution: " << std::endl
                  << " - middle box:     dx = " << dx_m << std::endl
                  << " - top/bottom box: dx = " << dx_t << std::endl;
  }
  if(lattice_type == 1) { // HCP lattice
    dy_m *= sqrt(3.) / 2.;
    dY_m = 2. * dy_m;
    dy_t *= sqrt(3.) / 2.;
    dY_t = 2. * dy_t;
    dz_m *= sqrt(2. / 3.);
    dZ_m = 2. * dz_m;
    dz_t *= sqrt(2. / 3.);
    dZ_t = 2. * dz_t;
    log_one(info) << "Lattice: HCP, resolution: " << std::endl
                  << " - middle box:     dx = " << dx_m << std::endl
                  << "                 2*dy = " << dY_m << std::endl
                  << "                 2*dz = " << dZ_m << std::endl
                  << " - top/bottom box: dx = " << dx_t << std::endl
                  << "                 2*dy = " << dY_t << std::endl
                  << "                 2*dz = " << dZ_t << std::endl;
  }
  if(lattice_type == 2) { // FCC lattice
    dy_m *= sqrt(3.) / 2.;
    dY_m = 2. * dy_m;
    dy_t *= sqrt(3.) / 2.;
    dY_t = 2. * dy_t;
    dz_m *= sqrt(2. / 3.);
    dZ_m = 3. * dz_m;
    dz_t *= sqrt(2. / 3.);
    dZ_t = 3. * dz_t;
    log_one(info) << "Lattice: FCC, resolution: " << std::endl
                  << " - middle box:     dx = " << dx_m << std::endl
                  << "                 2*dy = " << dY_m << std::endl
                  << "                 3*dz = " << dZ_m << std::endl
                  << " - top/bottom box: dx = " << dx_t << std::endl
                  << "                 2*dy = " << dY_t << std::endl
                  << "                 3*dz = " << dZ_t << std::endl;
  }

  // adjust width in y-direction of the middle block for symmetry
  w_m = floor(box_width / (3. * dY_m)) * dY_m;
  mbox_min[1] = -0.5 * w_m;
  mbox_max[1] = 0.5 * w_m + 0.01 * dy_m;

  // adjust top and bottom blocks
  gap = std::min(dY_m, dY_t) / 2.;
  if (lattice_type == 0) {
    gap *= 2.0;
  }
  w_t = floor((box_width / 2. - w_m / 2. - gap) / dY_t) * dY_t;
  tbox_min[1] = 0.5 * w_m + gap;
  tbox_max[1] = 0.5 * w_m + gap + w_t;
  bbox_min[1] = -0.5 * w_m - gap - w_t;
  bbox_max[1] = -0.5 * w_m - gap + 0.01 * dy_t;

  // set boxes length
  int Nx_t = floor(box_length / dx_t + 0.1);
  if(lattice_type == 0) {
    mbox_max[0] = box_length / 2.;
    mbox_min[0] = -box_length / 2. + dx_m / 2.;
    bbox_max[0] = Nx_t * dx_t / 2.;
    bbox_min[0] = -Nx_t * dx_t / 2. + dx_t / 2.;
    tbox_max[0] = Nx_t * dx_t / 2.;
    tbox_min[0] = -Nx_t * dx_t / 2. + dx_t / 2.;
  }
  else {
    mbox_max[0] = box_length / 2.;
    mbox_min[0] = -box_length / 2. + dx_m / 4.;
    bbox_max[0] = Nx_t * dx_t / 2.;
    bbox_min[0] = -Nx_t * dx_t / 2. + dx_t / 4.;
    tbox_max[0] = Nx_t * dx_t / 2.;
    tbox_min[0] = -Nx_t * dx_t / 2. + dx_t / 4.;
  }

  // set boxes height
  h_m = h_t = box_height;
  if constexpr(gdimension >= 3) {
    int Nz_m = floor(box_height / dZ_m);
    h_m = Nz_m * dZ_m;
    mbox_min[2] = -h_m / 2.;
    mbox_max[2] = h_m / 2.;

    int Nz_t = floor(box_height / dZ_t);
    h_t = Nz_t * dZ_t;
    bbox_min[2] = -h_t / 2.;
    bbox_max[2] = h_t / 2.;
    tbox_min[2] = -h_t / 2.;
    tbox_max[2] = h_t / 2.;
  }

  // count the number of particles
  np_middle = particle_lattice::count(
    lattice_type, domain_type, mbox_min, mbox_max, sph_separation, 0);
  np_top = particle_lattice::count(
    lattice_type, domain_type, tbox_min, tbox_max, sph_sep_t, np_middle);
  np_bottom = particle_lattice::count(lattice_type, domain_type, bbox_min,
    bbox_max, sph_sep_t, np_middle + np_top);

  SET_PARAM(nparticles, np_middle + np_bottom + np_top);
}

//----------------------------------------------------------------------------//
int
main(int argc, char * argv[]) {
  using namespace param;

  // launch MPI
  int rank, size;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  log_set_output_rank(0);

  // check options list: exactly one option is allowed
  if(argc != 2) {
    print_usage();
    MPI_Finalize();
    exit(0);
  }

  // only 2D and 3D cases are implemented
  assert(gdimension == 2 || gdimension == 3);

  // screen output
  log_one(info) << "Kelvin-Helmholtz instability initial data "
                << "in " << gdimension << "D" << std::endl;

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

  // erase ghosts                                                                                        
  auto pt = bodies.begin();                                                                              
  while (pt != bodies.end()) {                                                                           
    if (pt->type() == WALL) {                                                                            
      pt = bodies.erase(pt);                                                                             
    }                                                                                                    
    else {                                                                                               
      ++pt;                                                                                              
    }                                                                                                    
  }                                                                                                      
  if (modify_initial_data) {
    SET_PARAM(nparticles, bs.getNBodies());  
  }

  // screen output
  log_one(info) << "Number of particles: " << nparticles << std::endl;
  log_one(info) << "Initial data file: " << initial_data_file << std::endl;

  // allocate arrays
  // Position
  double * x = new double[nparticles]();
  double * y = new double[nparticles]();
  double * z = new double[nparticles]();

  if(not modify_initial_data) {
    // generate the lattice
    auto && [_npm, _npt, _npb] =
      std::make_tuple(particle_lattice::generate(lattice_type, domain_type,
                        mbox_min, mbox_max, sph_separation, 0, x, y, z),
        particle_lattice::generate(lattice_type, domain_type,
          tbox_min, tbox_max, sph_sep_t, np_middle, x, y, z),
        particle_lattice::generate(lattice_type, domain_type,
          bbox_min, bbox_max, sph_sep_t, nparticles - np_bottom, x, y, z));
    assert(np_middle == _npm && np_top == _npt && np_bottom == _npb);

    // stretch top and bottom blocks to align with the width
    double yx_stretch = floor(box_length / dx_t + 0.1) * dx_t / box_length;
    double yz_stretch = h_t / h_m;
    for(int a = np_middle; a < nparticles; ++a) {
      double y0 = w_m / 2. + gap;
      if(y[a] > 0) {
        y[a] = y0 + yx_stretch * yz_stretch * (y[a] - y0);
      }
      else {
        y[a] = -y0 + yx_stretch * yz_stretch * (y[a] + y0);
      }
      x[a] /= yx_stretch;
      if constexpr(gdimension == 3) {
        z[a] /= yz_stretch;
      }      
    }

    // stretch top and bottom blocks to align with the width
    double y_stretch = .5 * box_width / std::abs(y[np_middle + np_top]);
    double z_stretch = .5 * box_height / std::abs(z[np_middle + np_top]);
    for(int a = 0; a < nparticles; ++a) {
      y[a] *= y_stretch;
      if constexpr(gdimension == 3) {
        z[a] *= z_stretch;
      }
    }
    pmass *= y_stretch;
    if constexpr(gdimension == 3) {
      pmass *= z_stretch;
    }

    for(int64_t a = 0L; a < nparticles; ++a) {
      body & particle = bodies[a];
      if constexpr(gdimension == 1) {
        point_t pos = {x[a]};
        particle.set_coordinates(pos);
      }

      if constexpr(gdimension == 2) {
        point_t pos = {x[a], y[a]};
        particle.set_coordinates(pos);
      }

      if constexpr(gdimension == 3) {
        point_t pos = {x[a], y[a], z[a]};
        particle.set_coordinates(pos);
      }
    }      
  }

  // particle id number
  for(int64_t a = 0; a < nparticles; ++a) {
    body & particle = bodies[a];
    particle.set_id(a);
    double vy = 0.0;
    double vx = 0.0;      
    double vz = 0.0;
    point_t pos = particle.coordinates();
    if(std::abs(pos[1] - 0.25) < 0.025)
      vy = KH_A * sin(-2.0 * M_PI * (pos[0] + .5) / KH_lambda);
    if(std::abs(pos[1] + 0.25) < 0.025)
      vy = KH_A * sin( 2.0 * M_PI * (pos[0] + .5) / KH_lambda);
    if (std::abs(pos[1]) - 0.25 <= 0.0) {
      vx = vx_m;
    }
    else {
      vx = vx_t;
    }
    if(modify_initial_data) {
      if constexpr(gdimension == 1) {
        point_t vp = {vx};
        particle.setVelocity(vp);
      }
      if constexpr(gdimension == 2) {
        point_t vp = {vx,vy};
        particle.setVelocity(vp);
      }
      if constexpr(gdimension == 3) {    
        point_t vp = {vx,vy,vz};
        particle.setVelocity(vp);
      }
    }
    else {
      if(a < np_middle) {
        particle.setPressure(pressure_m);
        particle.setDensity(rho_m);
        particle.set_mass(pmass);
        if constexpr(gdimension == 1) {
          point_t vp = {vx_m};
          particle.setVelocity(vp);
        }
        if constexpr(gdimension == 2) {
          point_t vp = {vx_m,vy};
          particle.setVelocity(vp);
        }
        if constexpr(gdimension == 3) {
          point_t vp = {vx_m,vy,vz};
          particle.setVelocity(vp);
        }
        double u_a = pressure_m / (poly_gamma - 1.) / rho_m;
        particle.setInternalenergy(u_a);
        double h_a = sph_eta * kernels::kernel_width *
                     pow(pmass / rho_m, 1. / gdimension);
        particle.set_radius(h_a);
      }
      else {
        particle.setPressure(pressure_t);
        particle.setDensity(rho_t);
        particle.set_mass(pmass_t);
        if constexpr(gdimension == 1) {
          point_t vp = {vx_t};
          particle.setVelocity(vp);
        }
        if constexpr(gdimension == 2) {
          point_t vp = {vx_t,vy};
          particle.setVelocity(vp);
        }
        if constexpr(gdimension == 3) {
          point_t vp = {vx_t,vy,vz};
          particle.setVelocity(vp);
        }
        double u_a = pressure_t / (poly_gamma - 1.) / rho_t;
        particle.setInternalenergy(u_a);
        double h_a = sph_eta * kernels::kernel_width *
                     pow(pmass_t / rho_t, 1. / gdimension);
        particle.set_radius(h_a);
      }
    }
    particle.setDt(initial_dt);
  } // for part = 0 ... nparticles    

  log_one(info) << "Number of particles: " << nparticles << std::endl;

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

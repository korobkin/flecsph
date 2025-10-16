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
#include "influx.h"

#define SQ(x) ((x) * (x))
#define CU(x) ((x) * (x) * (x))
#define QU(x) ((x) * (x) * (x) * (x))

/*
Construct 3D compactified particle distribution for simulating realistic
kilonova ejecta, using an output flux through a spherical surface.
See 'influx.h' for details.
*/

//
// help message
//
void
print_usage() {
  std::cout << "Initial data generator for 3D kilonova ejecta from flux files\n" 
            << "or from an HDF5 file (such as the one made by Brendan)\n"
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

  // read input flux files
  influx::init();

  // reset spherical radius according to the median velocity and ejecta mass
  SET_PARAM(sphere_radius, influx::extraction_radius);
  log_one(info) << "flux extraction radius: " 
                << sphere_radius
                <<" [" << LENGTH_UNIT_STR << "]" << std::endl;

  // particle separation
  SET_PARAM(sph_separation, (2. * sphere_radius / (lattice_nx - 1)));

  total_mass = influx::total_ejecta_mass * phys::Msun; // convert to cgs
  rho_c = total_mass / CU(sphere_radius); // density estimate
  SET_PARAM(rho_initial, rho_c);
  log_one(info) << "average density: " << rho_c 
                << " [" << DENSITY_UNIT_STR << "]" << std::endl;

  // select equation of state and the type of lattice
  eos::select();

  // the value for initial timestep
  timestep = initial_dt;

  // bounding box of the domain
  bbox_min = -sphere_radius;
  bbox_max = sphere_radius;

  // total number of particles = Nx^3
  SET_PARAM(nparticles, CU(lattice_nx));

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
  sprintf(initial_data_file, "%s.h5part", initial_data_prefix);
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
 
  // set random seed (fix it for reproducibility)
  srand(12);
  //srand(time(0));

  // set simulation parameters
  param::mpi_read_params(argv[1]);
  set_derived_params();
  body_system<double, gdimension> bs;
  assert(not modify_initial_data); // not implemented
  bs.getLocalbodies().clear();
  bs.getLocalbodies().resize(nparticles);
  auto & bodies = bs.getLocalbodies();

  // Declare coordinate arrays
  double * x = new double[nparticles]();
  double * y = new double[nparticles]();
  double * z = new double[nparticles]();

  { using namespace influx;
    // Generate the particles layer-by-layer
    // 1. Create particle distribution over the time bins
    std::vector<size_t> Np_vs_time(INFLX_NT-1, 0);
    int Np_total = 0;
    double m1 = 0.0;
    for (int it=0; it<INFLX_NT-1; ++it) {
      double m2 = grid1d_cumulative_mass[it];
      Np_vs_time[it] = round((m2 - m1)*nparticles/total_ejecta_mass);
      Np_total += Np_vs_time[it];
      m1 = m2;
    }
    /* // randomly distribute over time slices
    int Np_total = nparticles;
    for (int it = 0; it < Np_total; ++it) {
      double x = (double)rand()/(double)RAND_MAX * total_ejecta_mass;
      auto j = interp::get_index(x, grid1d_cumulative_mass);
      Np_vs_time[j]++;
    }
    */

    // 2. Make the number of particles exact (it's not because of roundoff)
    int sgn = (Np_total < nparticles) ? 1 : -1;
    for (int64_t i = 0; i < std::abs((int64_t)nparticles 
                                   - (int64_t)Np_total); ++i) {
      double x = (double)rand()/(double)RAND_MAX * total_ejecta_mass;
      auto j = interp::get_index(x, grid1d_cumulative_mass);
      Np_vs_time[j] += sgn;
    }

    // 2,5. Compute vmax for every angular bin
    const double dphi = 2.*M_PI/INFLX_NPHI;
    int64_t a = 0L;
    // find the maximum velocity inside each solid angle bin
    double vr_max[INFLX_NTHETA][INFLX_NPHI];
    for (int ith=0; ith<INFLX_NTHETA; ith++) {
        for (int iph=0; iph<INFLX_NPHI; iph++) {
            double vr_tmp = 0.0;
            for (int it=0; it<INFLX_NT; it++) {
                vr_tmp = std::max(vr_tmp, grid3d_data[iph + INFLX_NPHI*(ith + INFLX_NTHETA*it)].vr);
            }
            vr_max[ith][iph] = vr_tmp;
        }
    }

    // 3. Distribute particles
    for (int it=1; it<INFLX_NT-1; ++it) {
      double * mass_it = grid3d_cumulative_mass.data() 
                       + it*INFLX_NTHETA*INFLX_NPHI;
      double m1 = grid1d_cumulative_mass[it-1];
      double m2 = grid1d_cumulative_mass[it];
      for (int i=0; i<Np_vs_time[it]; ++i) {
        grid_data_point_t gp;
        point_t pos;
        do {
          double x = m1 + (m2-m1)*(double)rand()/(double)RAND_MAX;
          auto ij = interp::get_index(x, mass_it, INFLX_NTHETA*INFLX_NPHI);

          int ith = ij / INFLX_NPHI;
          double theta = grid2d_theta[it*INFLX_NTHETA + ith];
          double dth = grid2d_dth[it*INFLX_NTHETA + ith];
          double c1 = cos(theta);
          double c2 = cos(theta + dth);
          c1 += (c2-c1)*(double)rand()/(double)RAND_MAX;
          double s1 = sqrt(1. - c1*c1);
          c1 = std::max(-1., std::min(1., c1));
          theta = acos(c1);
          
          int jphi = ij % INFLX_NPHI;
          double phi = (jphi + (double)rand()/(double)RAND_MAX)*dphi;
          
          pos = {s1*cos(phi), s1*sin(phi), c1};
          //auto gp = grid3d_data[jphi + INFLX_NPHI*(ith + INFLX_NTHETA*it)];
          double t = grid_times[it-1] + (grid_times[it] - grid_times[it-1])
                                        *(double)rand()/(double)RAND_MAX;
          gp = linear_interpolator(t, theta, phi);
          //hard code an expansion velocity 5x more than the largest velocity in flux files
          double rp = extraction_radius  + (5*vr_max[ith][jphi])*(grid_times[INFLX_NT] - t);
          pos *= rp;
          // calculate an adjusted density to account for large expansion
          gp.rho *= (extraction_radius*extraction_radius) / (rp*rp) * (gp.vr/(5*vr_max[ith][jphi]));
        } while (gp.vr <= 0. || std::isnan(gp.vr));
        bodies[a].set_id(a);
        bodies[a].set_coordinates(pos);
        bodies[a].set_state(INACTIVE);
        bodies[a].set_mass(mass_particle);
        bodies[a].setDensity(gp.rho);
        bodies[a].set_radius(cbrt(mass_particle/gp.rho));
        bodies[a].setAbar(initial_abar);
        bodies[a].setElectronfraction(initial_zbar/initial_abar);
        // bodies[a].setElectronfraction(gp.ye);
        bodies[a].setPressure(gp.pres);
        bodies[a].setTemperature(gp.temp);

        // Internal energy from nubhlight seems to be off; for now, 
        // compute internal energy from the ideal equation of state:
        // {rho, P} -> entropy -> internal energy
        //bodies[a].setInternalenergy(gp.uint);
        eos::compute_entropy(bodies[a]);
        eos::compute_internal_energy(bodies[a]);

        //point_t vel_r = {s1*cos(phi), s1*sin(phi), c1};
        //vel_r *= gp.vr;
        point_t vel = {gp.vr, gp.vth, gp.vphi};
        bodies[a].setVelocity(vel);
        ++a;
      }

    }

  } // using namespace influx

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

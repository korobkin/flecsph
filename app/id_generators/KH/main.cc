/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2017 Triad National Security, LLC
 * All rights reserved.
 *~--------------------------------------------------------------------------~*/


#include <iostream>
#include <algorithm>
#include <cassert>
#include <math.h>

#include "user.h"
#include "sodtube.h"
#include "params.h"
#include "lattice.h"
#include "kernels.h"
#include "io.h"

using namespace io;
//
// help message
//
void print_usage() {
  log_one(warn)
      << "Initial data generator for KH test in"
      << gdimension << "D" << std::endl
      << "Usage: ./KD_XD_generator <parameter-file.par>" << std::endl;
}

//
// derived parameters
//
static int64_t nparticlesproc;        // number of particles per proc
static double rho_m, rho_t;           // densities
static double vx_m, vx_t;             // velocities
static double pressure_m, pressure_t; // pressures
static std::string initial_data_file; // = initial_data_prefix + ".h5part"

// geometric extents of the three regions: top, middle and bottom
static point_t tbox_min, tbox_max;
static point_t mbox_min, mbox_max;
static point_t bbox_min, bbox_max;

// lattice spacing
// capital letters are for the full-periodicity lattice thickness
static double dx_m, dy_m, dY_m, dz_m, dZ_m;  // medium layer
static double dx_t, dy_t, dY_t, dz_t, dZ_t;  // medium layer

// width of the middle layer, top layer and gap between them
static double w_m, w_t, gap; 

static int64_t np_middle = 0; // number of particles in the middle block
static int64_t np_top = 0;    // number of particles in the top block
static int64_t np_bottom = 0; // number of particles in the bottom block
static double sph_sep_t = 0; // particle separation in top or bottom blocks
static double pmass = 0;      // particle mass in the middle block
static double pmass_t = 0;   // particle mass in top or bottom blocks

void set_derived_params() {
  using namespace std;
  using namespace param;

  // boundary tolerance factor
  const double b_tol = particle_lattice::b_tol;

  // support for only equal-mass configurations for now
  if(not equal_mass) {
    log_one(error) 
      << "Only equal-mass configurations are implemented"
      << std::endl;                    
    MPI_Finalize();
    exit(0);
  }

  mbox_max[0] = bbox_max[0] = tbox_max[0] = box_length/2.;
  mbox_min[0] = bbox_min[0] = tbox_min[0] =-box_length/2.;

  mbox_max[1] =  box_width/4.;
  mbox_min[1] = -box_width/4.;
  bbox_min[1] = -box_width/2.;
  bbox_max[1] = -box_width/4.;
  tbox_min[1] =  box_width/4.;
  tbox_max[1] =  box_width/2.;

  if(gdimension == 3){
    mbox_max[2] = bbox_max[2] = tbox_max[2] = box_height/2.;
    mbox_min[2] = bbox_min[2] = tbox_min[2] =-box_height/2.;
  }


  // set physical parameters
  // --  in the top and bottom boxes (tbox and bbox):
  rho_t      =  rho_initial;
  pressure_t =  pressure_initial;
  vx_t       = -flow_velocity/2.0;
  // -- in the middle box
  rho_m      =  rho_t * density_ratio; // 2.0 by default
  pressure_m =  pressure_t;          // pressures must be equal in KH test
  vx_m       =  flow_velocity/2.0;

  // file to be generated
  std::ostringstream oss;
  oss << initial_data_prefix << ".h5part";
  initial_data_file = oss.str();

  // select particle lattice and kernel function
  particle_lattice::select();
  kernels::select();

  // particle mass and spacing
  SET_PARAM(sph_separation, box_length/lattice_nx);
  if(gdimension == 3){
    pmass = rho_m*sph_separation*sph_separation*sph_separation;
    if (lattice_type == 1 or lattice_type == 2)
      pmass *= 1./sqrt(2.);
    sph_sep_t = sph_separation*cbrt(density_ratio);
  }
  if(gdimension == 2){
    pmass = rho_m*sph_separation*sph_separation;
    if (lattice_type == 1 or lattice_type == 2)
      pmass *= sqrt(3)/2;
    sph_sep_t = sph_separation*sqrt(density_ratio);
  }
  pmass_t = pmass;

  // lattice spacing
  dx_m = dy_m = dY_m = dz_m = dZ_m = sph_separation;
  dx_t = dy_t = dY_t = dz_t = dZ_t = sph_sep_t;
  if (lattice_type == 0) {
    log_one(info)
      << "Lattice: rectangular, resolution: " << std::endl
      << " - middle box:     dx = " << dx_m << std::endl
      << " - top/bottom box: dx = " << dx_t << std::endl;
  }
  if (lattice_type == 1) { // HCP lattice
    dy_m *= sqrt(3.)/2.;  dY_m = 2.*dy_m;       
    dy_t *= sqrt(3.)/2.;  dY_t = 2.*dy_t;
    dz_m *= sqrt(2./3.);  dZ_m = 2.*dz_m; 
    dz_t *= sqrt(2./3.);  dZ_t = 2.*dz_t; 
    log_one(info)
      << "Lattice: HCP, resolution: " << std::endl
      << " - middle box:     dx = " << dx_m << std::endl
      << "                 2*dy = " << dY_m << std::endl
      << "                 2*dz = " << dZ_m << std::endl
      << " - top/bottom box: dx = " << dx_t << std::endl
      << "                 2*dy = " << dY_t << std::endl
      << "                 2*dz = " << dZ_t << std::endl;
  }
  if (lattice_type == 2) { // FCC lattice
    dy_m *= sqrt(3.)/2.;  dY_m = 2.*dy_m;       
    dy_t *= sqrt(3.)/2.;  dY_t = 2.*dy_t;
    dz_m *= sqrt(2./3.);  dZ_m = 3.*dz_m; 
    dz_t *= sqrt(2./3.);  dZ_t = 3.*dz_t; 
    log_one(info)
      << "Lattice: FCC, resolution: " << std::endl
      << " - middle box:     dx = " << dx_m << std::endl
      << "                 2*dy = " << dY_m << std::endl
      << "                 3*dz = " << dZ_m << std::endl
      << " - top/bottom box: dx = " << dx_t << std::endl
      << "                 2*dy = " << dY_t << std::endl
      << "                 3*dz = " << dZ_t << std::endl;
  }

  // adjust width in y-direction of the middle block for symmetry
  w_m = floor(box_width/(3.*dY_m)) * dY_m;
  mbox_min[1] = -0.5*w_m;
  mbox_max[1] =  0.5*w_m + 0.5*dy_m;

  // adjust top and bottom blocks
  gap = (dy_m + dy_t)/2.; //std::min(dY_m,dY_t)/2.;
  w_t = floor((box_width/2. - w_m/2. - gap)/dY_t) * dY_t;
log_one(info) << std::endl 
  << "w_m = " << w_m << std::endl
  << "w_t = " << w_t << std::endl
  << "gap = " << gap << std::endl
  << "box_width before = " << box_width << std::endl;

  //SET_PARAM(box_width, w_m + 2.*(w_t + gap))
log_one(info) << std::endl 
  << "box_width after  = " << box_width << std::endl
  << "w_m+2.*(w_t+gap) = " << w_m+2.*(w_t+gap) << std::endl;
  tbox_min[1] = 0.5*w_m + gap;
  tbox_max[1] = 0.5*w_m + gap + w_t;
  bbox_min[1] =-0.5*w_m - gap - w_t;
  bbox_max[1] =-0.5*w_m - gap + 0.5*dy_t;

  // set boxes length 
  int Nx_t = floor(box_length/dx_t + 0.1);
  if (lattice_type == 0) {
    mbox_max[0] =  box_length/2.;
    mbox_min[0] = -box_length/2. + dx_m/2.;
    bbox_max[0] =  Nx_t*dx_t/2.;
    bbox_min[0] = -Nx_t*dx_t/2. + dx_t/2.;
    tbox_max[0] =  Nx_t*dx_t/2.;
    tbox_min[0] = -Nx_t*dx_t/2. + dx_t/2.;
  }
  else {
    mbox_max[0] =  box_length/2.;
    mbox_min[0] = -box_length/2. + dx_m/4.;
    bbox_max[0] =  Nx_t*dx_t/2.;
    bbox_min[0] = -Nx_t*dx_t/2. + dx_t/4.;
    tbox_max[0] =  Nx_t*dx_t/2.;
    tbox_min[0] = -Nx_t*dx_t/2. + dx_t/4.;
  }

  // set boxes height (TODO)
  if constexpr (gdimension >= 3) {
    int Nz = (int)(box_height/dZ_m) - 1;
    for(int k = Nz; k < Nz*100; ++k) {
      double w2 = floor(k*dZ_m/dZ_t + b_tol)*dZ_t;
      if (fabs(w2 - k*dZ_m) < lattice_mismatch_tolerance*dZ_t) {
        SET_PARAM(box_height, std::min(w2,k*dZ_m));
        bbox_min[2] = -box_height/2.;
        bbox_max[2] =  box_height/2.;
        mbox_min[2] = -box_height/2.;
        mbox_max[2] =  box_height/2.;
        tbox_min[2] = -box_height/2.;
        tbox_max[2] =  box_height/2.;
        break;
      }
    }
  }

  // report adjusted dimensions
  log_one(warn) 
    << "Domain has been adjusted for periodic boundaries." << std::endl;
  log_one(warn) 
    << "For evolution, modify domain dimensions as follows:" << std::endl 
    << "  box_length = " << box_length << std::endl 
    << "  box_width = "  << box_width  << std::endl 
    << "  box_height = " << box_height << std::endl;


  log_one(warn) << "Lattice mismatch, X-direction:" << std::endl
      << " -    top box: "
      <<   (box_length - floor((tbox_max[0]-tbox_min[0])/dx_t)*dx_t) 
      <<   ", dx = " << dx_t << ", mismatch/dx = " 
      <<   (box_length/dx_t - floor((tbox_max[0]-tbox_min[0])/dx_t)) 
      << std::endl
      << " - middle box: "
      <<   (box_length - floor((mbox_max[0]-mbox_min[0])/dx_m)*dx_m) 
      <<   ", dx = " << dx_m << ", mismatch/dx = " 
      <<   (box_length/dx_m - floor((mbox_max[0]-mbox_min[0])/dx_m)) 
      << std::endl
      << " - bottom box: "
      <<   (box_length - floor((bbox_max[0]-bbox_min[0])/dx_t)*dx_t) 
      <<   ", dx = " << dx_t << ", mismatch/dx = " 
      <<   (box_length/dx_t - floor((bbox_max[0]-bbox_min[0])/dx_t)) 
      << std::endl;

  if constexpr (gdimension >= 3) 
    log_one(warn) << "Lattice mismatch, Z-direction:" << std::endl
      << " -    top box: "
      <<   (box_height - floor((tbox_max[2]-tbox_min[2])/dZ_t)*dZ_t) 
      <<   ", dz = " << dZ_t << ", mismatch/dz = " 
      <<   (box_height/dZ_t - floor((bbox_max[2]-bbox_min[2])/dZ_t)) 
      << std::endl
      << " - middle box: "
      <<   (box_height - floor((mbox_max[2]-mbox_min[2])/dZ_m)*dZ_m) 
      <<   ", dz = " << dZ_m << ", mismatch/dz = " 
      <<   (box_height/dZ_m - floor((mbox_max[2]-mbox_min[2])/dZ_m)) 
      << std::endl
      << " - bottom box: "
      <<   (box_height - floor((bbox_max[2]-bbox_min[2])/dZ_t)*dZ_t) 
      <<   ", dz = " << dZ_t << ", mismatch/dz = " 
      <<   (box_height/dZ_t - floor((bbox_max[2]-bbox_min[2])/dZ_t)) 
      << std::endl;

  // count the number of particles
  np_middle = particle_lattice::count(lattice_type,gdimension,mbox_min,mbox_max,
                                      sph_separation,0);
  np_top    = particle_lattice::count(lattice_type,gdimension,tbox_min,tbox_max,
                                      sph_sep_t, np_middle);
  np_bottom = particle_lattice::count(lattice_type,gdimension,bbox_min,bbox_max,
                                      sph_sep_t, np_middle + np_top);

  SET_PARAM(nparticles, np_middle + np_bottom + np_top);

}

//----------------------------------------------------------------------------//
int main(int argc, char * argv[]){
  using namespace param;

  // launch MPI
  int rank, size;
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  log_set_output_rank(0);

  // check options list: exactly one option is allowed
  if (argc != 2) {
    print_usage();
    MPI_Finalize();
    exit(0);
  }

  // anything other than 2D is not implemented yet
  assert (gdimension == 2 || gdimension == 3);
  assert (domain_type == 0);

  // screen output
  log_one(info)
    << "Kelvin-Helmholtz instability initial data " 
    << "in " << gdimension << "D" << std::endl;

  // set simulation parameters
  param::mpi_read_params(argv[1]);
  set_derived_params();

  // screen output
  log_one(info)
    << "Number of particles: "
    << nparticles << std::endl;
  log_one(info)
    << "Initial data file: " 
    << initial_data_file << std::endl;

  // allocate arrays
  // Position
  double* x = new double[nparticles]();
  double* y = new double[nparticles]();
  double* z = new double[nparticles]();
  // Velocity
  double* vx = new double[nparticles]();
  double* vy = new double[nparticles]();
  double* vz = new double[nparticles]();
  // Acceleration
  double* ax = new double[nparticles]();
  double* ay = new double[nparticles]();
  double* az = new double[nparticles]();
  // Smoothing length
  double* h = new double[nparticles]();
  // Density
  double* rho = new double[nparticles]();
  // Internal Energy
  double* u = new double[nparticles]();
  // Pressure
  double* P = new double[nparticles]();
  // Mass
  double* m = new double[nparticles]();
  // Id
  int64_t* id = new int64_t[nparticles]();
  // Timestep
  double* dt = new double[nparticles]();

  // generate the lattice
  assert (np_middle == particle_lattice::generate( lattice_type,gdimension,
          mbox_min,mbox_max,sph_separation,0,x,y,z));
  assert (np_top    == particle_lattice::generate( lattice_type,gdimension,
          tbox_min,tbox_max,sph_sep_t,np_middle,x,y,z));
  assert (np_bottom == particle_lattice::generate( lattice_type,gdimension,
          bbox_min,bbox_max,sph_sep_t,nparticles-np_bottom,x,y,z));

  // stretch top and bottom blocks to align with the width
  double yx_stratch = floor(box_length/sph_sep_t + 0.1)*sph_sep_t/box_length;
  double dx_plus = 0.5*(box_length - yx_stratch*box_length);
log_one(info) << "yx_stratch = " << yx_stratch << std::endl;
log_one(info) << "dx_plus =    " << dx_plus << std::endl;
  for (int i = np_middle; i < nparticles; ++i) {
    y[i] *= yx_stratch;
    x[i] /= yx_stratch;
  }
    
  // stretch top and bottom blocks to align with the width
  double y_stratch = .5*box_width / std::abs(y[np_middle + np_top]);
  for (int i = 0; i < nparticles; ++i) {
    y[i] *= y_stratch;
  }
  pmass *= y_stratch;

  // max. value for the speed of sound
  double cs = sqrt(poly_gamma*std::max(pressure_m/rho_m,pressure_t/rho_t));

  // suggested timestep
  double timestep = timestep_cfl_factor
                  * sph_separation/std::max(cs, flow_velocity);

  // particle id number
  int64_t posid = 0;
  double wmid2 = -y[0] + .1*dy_m;
  for(int64_t part=0; part<nparticles; ++part){
    id[part] = posid++;
    if (std::abs(y[part]) < wmid2) {
      P[part] = pressure_m;
      rho[part] = rho_m;
      vx[part] = vx_m;
      m[part] = pmass;
    }
    else {
      P[part] = pressure_t;
      rho[part] = rho_t;
      vx[part] = vx_t;
      m[part] = pmass_t;
    }

    vy[part] = 0.;

    // Add velocity perturbation a-la Price (2008)
    if(y[part] < 0.25 and y[part] > 0.25 - 0.025)
      vy[part] = KH_A*sin(-2*M_PI*(x[part]+.5)/KH_lambda);
    if(y[part] >-0.25 and y[part] <-0.25 + 0.025)
      vy[part] = KH_A*sin(2*M_PI*(x[part]+.5)/KH_lambda);

    // compute internal energy using gamma-law eos
    u[part] = P[part]/(poly_gamma-1.)/rho[part];

    // particle masses and smoothing length
    m[part] = pmass;
    h[part] = sph_eta * kernels::kernel_width
                      * pow(m[part]/rho[part],1./gdimension);

  } // for part=0..nparticles

  // delete the output file if exists
  remove(initial_data_file.c_str());

  hid_t dataFile = H5P_openFile(initial_data_file.c_str(),H5F_ACC_RDWR);

  int use_fixed_timestep = 1;
  // add the global attributes
  H5P_writeAttribute(dataFile,"nparticles",&nparticles);
  H5P_writeAttribute(dataFile,"timestep",&timestep);
  int dim = gdimension;
  H5P_writeAttribute(dataFile,"dimension",&dim);
  H5P_writeAttribute(dataFile,"use_fixed_timestep",&use_fixed_timestep);

  H5P_setNumParticles(nparticles);
  H5P_setStep(dataFile,0);

  //H5PartSetNumParticles(dataFile,nparticles);
  H5P_writeDataset(dataFile,"x",x,nparticles);
  H5P_writeDataset(dataFile,"y",y,nparticles);
  H5P_writeDataset(dataFile,"z",z,nparticles);
  H5P_writeDataset(dataFile,"vx",vx,nparticles);
  H5P_writeDataset(dataFile,"vy",vy,nparticles);
  H5P_writeDataset(dataFile,"h",h,nparticles);
  H5P_writeDataset(dataFile,"rho",rho,nparticles);
  H5P_writeDataset(dataFile,"u",u,nparticles);
  H5P_writeDataset(dataFile,"P",P,nparticles);
  H5P_writeDataset(dataFile,"m",m,nparticles);
  H5P_writeDataset(dataFile,"id",id,nparticles);

  H5P_closeFile(dataFile);

  delete[] x, y, z, vx, vy, vz, ax, ay, az, h, rho, u, P, m, id, dt;

  MPI_Finalize();
  return 0;
}

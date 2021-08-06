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

#include <hdf5.h>
#define IND3D(i,j,k) ((i) + (j)*(Nx + (k)*Nz))
/*
  Regridding: reading an input h5part file with particles and interpolating
  it on a grid, using the current kernel and parameters.
*/

//
// help message
//
void
print_usage() {
  std::cout << "Gridding particle data\n"
            << "Usage: ./sedov_generator <parameter-file.par>\n";
}

//
// derived parameters
//
static double timestep = 1.0;       // Recommended timestep
static double total_mass = 1.;      // total mass of the fluid
static double mass_particle = 1.;   // mass of an individual particle
static point_t bbox_max, bbox_min;  // bounding box of the domain
static char initial_data_file[256]; // = initial_data_prefix[_XXXXX].h5part"
static double dx = 0.1;             // grid step
static int Nx = 150;
static int Ny = 150;
static int Nz = 150;

void
set_derived_params() {
  using namespace param;

  eos::select();
  density_profiles::select();
  particle_lattice::select();

  // The value for constant timestep
  timestep = initial_dt;

  // Bounding box
  bbox_min[0] = -box_length / 2.;
  bbox_max[0] =  box_length / 2.;
  bbox_min[1] = -box_width / 2.;
  bbox_max[1] =  box_width / 2.;
  bbox_min[2] = -box_height / 2.;
  bbox_max[2] =  box_height / 2.;

  // particle separation
  Nx = lattice_nx;
  dx = box_length / (double)(lattice_nx - 1);
  Ny = (int)(box_width / dx) + 1;
  Nz = (int)(box_height / dx) + 1;

  // set kernel
  kernels::select();

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
  double pt_mass = bodies[0].mass();

  // Declare coordinate arrays
  double * x = new double[nparticles]();
  double * y = new double[nparticles]();
  double * z = new double[nparticles]();

  // Declare grid coordinates
  double *xg = new double[Nx]();
  double *yg = new double[Ny]();
  double *zg = new double[Nz]();

  // Rescale everything to this time
  double zoom_factor = 86400.0/6.510414752339e+01;

  // Create the coordinate data.
  for(int i=0; i<Nx; ++i) xg[i] = -.5*box_length + dx*i;
  for(int j=0; j<Ny; ++j) yg[j] = -.5*box_width  + dx*j;
  for(int k=0; k<Nz; ++k) zg[k] = -.5*box_height + dx*k;

  // simple density estimation
  int ndx = 0;
  const int64_t Nxyz = Nx*Ny*Nz;
  double * rho = new double[Nxyz]();
  double * ye  = new double[Nxyz]();
  memset(rho, 0x00, sizeof(double)*Nxyz);
  memset(ye,  0x00, sizeof(double)*Nxyz);
  point_t bbox_min{0}, bbox_max{0};
  for(int64_t a = 0L; a < nparticles; ++a) {
    const point_t rp = bodies[a].coordinates();
    for (int d=0; d<3; ++d) {
      if (rp[d] < bbox_min[d]) bbox_min[d] = rp[d];
      if (rp[d] > bbox_max[d]) bbox_max[d] = rp[d];
    }

    int ix = floor((rp[0] - xg[0])/dx);
    if (ix > Nx - 1 || ix < 0) continue;

    int jy = floor((rp[1] - yg[0])/dx);
    if (jy > Ny - 1 || jy < 0) continue;

    int kz = floor((rp[2] - zg[0])/dx);
    if (kz > Nz - 1 || kz < 0) continue;

    int ijk = ix + Nx*(jy + Ny*kz);
    rho[ijk] += pt_mass / (dx*dx*dx);
    ye[ijk] = bodies[a].getElectronfraction();
  }
  std::cout << bbox_min << " : " << bbox_max << "\n";

  // a better interpolator
  // 1. determine the smoothing length for every grid point
  double * hg  = new double[Nxyz]();  // smoothing length for the grid points
  double * dgp = new double[Nxyz]();  // distance to the nearest particle
  int64_t * np = new int64_t[Nxyz]();  // index of the nearest particle
  memset(hg, 0x00, sizeof(double)*Nxyz);
  memset(dgp, 0x00,sizeof(double)*Nxyz);
  memset(np,0x00, sizeof(int64_t)*Nxyz);
  for(int64_t a = 0L; a < nparticles; ++a) {
    double h_a = bodies[a].radius();
    point_t rp = bodies[a].coordinates();
    int imn = floor((rp[0] - h_a - xg[0])/dx);
    if (imn > Nx - 1) continue;
    imn = std::max(imn, 0);
    int imx = floor((rp[0] + h_a - xg[0])/dx);
    if (imx < 0) continue;
    imx = std::min(imx, Nx-1);

    int jmn = floor((rp[1] - h_a - yg[0])/dx);
    if (jmn > Ny - 1) continue;
    jmn = std::max(jmn, 0);
    int jmx = floor((rp[1] + h_a - yg[0])/dx);
    if (jmx < 0) continue;
    jmx = std::min(jmx, Ny-1);

    int kmn = floor((rp[2] - h_a - zg[0])/dx);
    if (kmn > Nz - 1) continue;
    kmn = std::max(kmn, 0);
    int kmx = floor((rp[2] + h_a - zg[0])/dx);
    if (kmx < 0) continue;
    kmx = std::min(kmx, Nz-1);

    for(int i=imn;i<imx;++i)
    for(int j=jmn;j<jmx;++j)
    for(int k=kmn;k<kmx;++k) {
      double r = sqrt(SQ(xg[i]-rp[0]) + SQ(yg[j]-rp[1]) + SQ(zg[k]-rp[2]));
      if (r > h_a) continue;
      int64_t ijk = i + Nx*(j + Ny*k);
      if (np[ijk] == 0     // first time for this grid point
      ||  r < dgp[ijk]) {  // or the particle is closer than previous ones
        hg[ijk] = h_a;
        dgp[ijk] = r;
        np[ijk] = a;
        continue;
      }
    }
  } // for a...

  // just take density and electron fraction from the neares point
  for(int i=0;i<Nx;++i)
  for(int j=0;j<Ny;++j)
  for(int k=0;k<Nz;++k) {
    int64_t ijk = i + Nx*(j + Ny*k);
    int64_t a = np[ijk];
    if (a == 0) continue;
    body & pt = bodies[a];
    rho[ijk] = pt.getDensity()/CU(zoom_factor);
    ye[ijk] = pt.getElectronfraction();
  }

  // Rescale coordinates
  for(int i=0; i<Nx; ++i) xg[i] *= zoom_factor;
  for(int j=0; j<Ny; ++j) yg[j] *= zoom_factor;
  for(int k=0; k<Nz; ++k) zg[k] *= zoom_factor;

  // open the hdf5 file
  hid_t     file_id;
  char h5fname[1024];
  sprintf(h5fname, "%s.h5", output_h5data_prefix);
  file_id = H5Fcreate(h5fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  // create dataspace / dataset
  hid_t     dataset_id, dataspace_id, status;
  hsize_t   dims[3];
  dims[0] = Nx;
  dims[1] = Ny;
  dims[2] = Nz;

  // record x-grid
  dataspace_id = H5Screate_simple(1, &(dims[0]), NULL);
  dataset_id = H5Dcreate(file_id, "/x", H5T_NATIVE_DOUBLE, dataspace_id,
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                    H5P_DEFAULT, xg);
  status = H5Dclose(dataset_id);
  status = H5Sclose(dataspace_id);

  // record y-grid
  dataspace_id = H5Screate_simple(1, &(dims[1]), NULL);
  dataset_id = H5Dcreate(file_id, "/y", H5T_NATIVE_DOUBLE, dataspace_id,
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                    H5P_DEFAULT, yg);
  status = H5Dclose(dataset_id);
  status = H5Sclose(dataspace_id);

  // record z-grid
  dataspace_id = H5Screate_simple(1, &(dims[2]), NULL);
  dataset_id = H5Dcreate(file_id, "/z", H5T_NATIVE_DOUBLE, dataspace_id,
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                    H5P_DEFAULT, zg);
  status = H5Dclose(dataset_id);
  status = H5Sclose(dataspace_id);

  // record the density
  dataspace_id = H5Screate_simple(3, dims, NULL);
  dataset_id = H5Dcreate(file_id, "/rho", H5T_NATIVE_DOUBLE,
                         dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,
                    H5S_ALL, H5P_DEFAULT, rho);
  status = H5Dclose(dataset_id);
  status = H5Sclose(dataspace_id);

  // electron fraction
  dataspace_id = H5Screate_simple(3, dims, NULL);
  dataset_id = H5Dcreate(file_id, "/ye", H5T_NATIVE_DOUBLE,
                         dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,
                    H5S_ALL, H5P_DEFAULT, ye);
  status = H5Dclose(dataset_id);
  status = H5Sclose(dataspace_id);

  // xdmf file
  char xdmf_out[1024];
  sprintf(xdmf_out, "%s.xdmf", output_h5data_prefix);
  remove(xdmf_out);

  std::ostringstream xdmf_text;
  xdmf_text <<
      "<?xml version=\"1.0\" ?>\n"
      "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\">\n"
      "<Xdmf Version=\"2.0\">\n"
      " <Domain>\n"
      "   <Grid Name=\"sph-regridded\" GridType=\"Uniform\">\n"
      "     <!-- dimensions: NZ NY NX -->\n"
      "     <Topology TopologyType=\"3DRectMesh\" "
              "Dimensions=\""<< Nz <<" "<< Ny <<" "<< Nx <<"\"/>\n"
      "     <Geometry GeometryType=\"VXVYVZ\">\n"
      "       <DataItem Format=\"HDF\" Dimensions=\""<< Nx
              <<"\" NumberType=\"Float\" Precision=\"4\">\n"
      "           "<< output_h5data_prefix <<".h5:/x\n"
      "       </DataItem>\n"
      "       <DataItem Format=\"HDF\" Dimensions=\""<< Ny
              <<"\" NumberType=\"Float\" Precision=\"4\">\n"
      "           "<< output_h5data_prefix <<".h5:/y\n"
      "       </DataItem>\n"
      "       <DataItem Format=\"HDF\" Dimensions=\""<< Nz
              <<"\" NumberType=\"Float\" Precision=\"4\">\n"
      "           "<< output_h5data_prefix <<".h5:/z\n"
      "       </DataItem>\n"
      "     </Geometry>\n"
      "     <Attribute Name=\"rho\" Center=\"Node\">\n"
      "       <!-- dimensions: NX NY NZ -->\n"
      "       <DataItem Format=\"HDF\" "
                  "Dimensions=\""<< Nx <<" "<< Ny <<" "<< Nz
                  <<"\" NumberType=\"Float\" Precision=\"4\">\n"
      "           "<< output_h5data_prefix <<".h5:/rho\n"
      "       </DataItem>\n"
      "     </Attribute>\n"
      "     <Attribute Name=\"ye\" Center=\"Node\">\n"
      "       <DataItem Format=\"HDF\" "
                  "Dimensions=\""<< Nx <<" "<< Ny <<" "<< Nz
                  <<"\" NumberType=\"Float\" Precision=\"4\">\n"
      "           "<< output_h5data_prefix <<".h5:/ye\n"
      "       </DataItem>\n"
      "     </Attribute>\n"
      "   </Grid>\n"
      " </Domain>\n"
      "</Xdmf>\n";

  // Open file in append mode
  remove(xdmf_out);
  std::ofstream out(xdmf_out, std::ios_base::app);
  out << xdmf_text.str();
  out.close();

  //bs.write_bodies(initial_data_prefix, 0, 0.0);
  MPI_Finalize();

  // cleanup
  delete[] x, y, z, xg, yg, zg, rho, ye;

  return 0;
}

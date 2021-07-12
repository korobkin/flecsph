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
 * @brief Specialization and Main driver used in FleCSI.
 * The Specialization Driver is normally used to register data and the main
 * code is in the Driver.
 */

#include <iostream>
#include <numeric> // For accumulate

#include <mpi.h>
#include <omp.h>

// #define poly_gamma 5./3.
#include "analysis.h"
#include "bodies_system.h"
#include "default_physics.h"
#include "diagnostic.h"
#include "params.h"

#include "control.h"
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

  // iteration and time
  physics::iteration = initial_iteration;
  physics::totaltime = initial_time;
  physics::dt = initial_dt;

  // set equation of state
  eos::select();

  // set external force
  external_force::select(external_force_type);
}

int
advance() {
  using namespace param;

  auto& parameter_file = control::policy().filename();


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

  size_t total = 5000;
  do {
    log_one(info)<<"######## Iteration: "<<total<<std::endl; 
    MPI_Barrier(MPI_COMM_WORLD); 
    //analysis::screen_output(rank);
    bs.update_iteration();
    double begin = omp_get_wtime();
    size_t total = 0;
    
    bs.apply_in_smoothinglength(
      [&](tree_topology_t::entity_t & e,
        std::vector<tree_topology_t::entity_t *> & n, size_t & total) {
        total += n.size();
        bool found = false;
        for(auto nb : n) {
          e.id() == nb->id() ? found = true : found;
        }
        assert(found);
      },
      total);
    std::cout << "Average: " << total / bs.nbodies() << std::endl;
    double end = omp_get_wtime();
    std::cout << "Traversal time: " << end - begin << "s " << std::endl;

#if 0 
    bs.reset_ghosts(); 
    begin = omp_get_wtime(); 
    total = 0; 
    bs.apply_in_smoothinglength(
      [&](tree_topology_t::entity_t& e, std::vector<tree_topology_t::entity_t*> & n, size_t& total){
        total+=n.size();
        bool found = false; 
        auto id_e = e.id(); 
        for(auto nb: n){
          e.id() == nb->id()?found=true:found; 
        } 
        assert(found); 
      },total
    );
    std::cout<<"Average 2: "<<total/bs.nbodies()<<std::endl;
    end = omp_get_wtime(); 
    std::cout<<"Traversal time 2: "<<end-begin<<"s "<<std::endl;
#endif
    ++physics::iteration;
  } while(--total != 0);
  return 0; 
} // advance

bool
check_conservation(const std::vector<analysis::e_conservation> & check) {
  return analysis::check_conservation(check);
}

control::action<advance, cp::advance> advance_action;

int
main(int argc, char * argv[]) {

  auto status = flecsi::initialize(argc, argv);
  auto pf = parameter_file.value(); 
  status = control::check_status(status);
  if(status != flecsi::run::status::success) {
    return status < flecsi::run::status::clean ? 0 : status;
  }
  flecsi::log::add_output_stream("clog", std::clog, true);

  auto& filename = control::policy().filename();
  filename = pf; 

  status = flecsi::start(control::execute);
  flecsi::finalize();
  return status;
}

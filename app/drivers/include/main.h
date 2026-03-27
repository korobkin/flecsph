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

#pragma once

#include <iostream>
#include <string>
#include <mpi.h>

// Forward declaration: each driver defines this function
int advance(const std::string& parameter_file);

inline int
driver_main(int argc, char * argv[]) {
  MPI_Init(&argc, &argv);
  if(argc < 2) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if(rank == 0) {
      std::cerr << "Usage: " << argv[0] << " <parameters.par>" << std::endl;
    }
    MPI_Finalize();
    return 1;
  }
  std::string parameter_file(argv[1]);
  if(parameter_file.find(".par") == std::string::npos) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if(rank == 0) {
      std::cerr << "Error: file(" << parameter_file
                << ") has invalid suffix (expected .par)" << std::endl;
    }
    MPI_Finalize();
    return 1;
  }
  int status = advance(parameter_file);
  MPI_Finalize();
  return status;
}

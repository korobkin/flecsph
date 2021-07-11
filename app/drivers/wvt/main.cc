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
 * @file main.cc
 * @author Julien Loiseau
 * @date April 2017
 * @brief Main function, start MPI with Gasnet. Then launch fleCSI runtime.
 */

#include "flecsi/execution.hh"
#include "control.h"

void
usage(char progname[]) {
  std::cout << "Usage: " << progname << " <parameter-file.par>" << std::endl;
}

int
main(int argc, char * argv[]) {

  /*
    Check options list for a positional option which represents parameter file
    remove it to pass the remaining options to FleCSI
   */
  char * parameter_file;
  int i = 1;
  for (; i < argc; ++i) {
    if (argv[i][0] != '-') {
      parameter_file = (char*)(argv[i]);
      for (int j = i; j < argc - 1; ++j) {
        argv[j] = argv[j+1];
      }
      break;
    }
  }
  if(i == argc) {
    std::cerr << "ERROR: parameter file not specified!" << std::endl;
    usage(argv[0]);
    return 1;
  }

  --argc; // we removed one argument
  auto status = flecsi::initialize(argc, argv);
  status = control::check_status(status);
  if(status != flecsi::run::status::success) {
    return status < flecsi::run::status::clean ? 0 : status;
  }
  flecsi::log::add_output_stream("clog", std::clog, true);

  auto& filename = control::policy().filename();
  filename = parameter_file; 

  /*
    Pass the control model's 'execute' method to start. FleCSI will invoke
    the execute function after runtime initialization. This will, in turn,
    execute all of the cycles, and actions of the control model.
   */
  status = flecsi::start(control::execute);
  flecsi::finalize();
  return status;
}

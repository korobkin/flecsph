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
auto status = flecsi::initialize(1, argv);
  /*
    The check_options() method checks to see if any control-model options were
    specified on the command line, and handles them appropriately.
   */

  status = control::check_status(status);
  if(status != flecsi::run::status::success) {
    return status < flecsi::run::status::clean ? 0 : status;
  }
  flecsi::log::add_output_stream("clog", std::clog, true);

  // check options list: exactly one option is allowed
  if(argc != 2) {
    std::cerr << "ERROR: parameter file not specified!" << std::endl;
    usage(argv[0]);
    return 1;
  }
  auto& filename = control::policy().filename();
  filename = argv[1]; 

  /*
    Pass the control model's 'execute' method to start. FleCSI will invoke
    the execute function after runtime initialization. This will, in turn,
    execute all of the cycles, and actions of the control model.
   */
  status = flecsi::start(control::execute);
  flecsi::finalize();
  return status;
}

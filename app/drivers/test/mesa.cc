#include "gtest/gtest.h"

#include <cmath>
#include <iostream>
#include <flecsi/flog.hh>
#include "control.h"

#include <mpi.h>

using namespace flecsi;

int
main(int argc, char * argv[]) {
  auto status = flecsi::initialize(argc, argv);
  status = control::check_status(status);
  if(status != flecsi::run::status::success) {
    return status < flecsi::run::status::clean ? 0 : status;
  }
  control::policy().filename() = "mesa_nx20.par";
  flecsi::log::add_output_stream("clog", std::clog, true);
  status = flecsi::start(control::execute);
  flecsi::finalize();
  return status;
}

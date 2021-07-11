#include "gtest/gtest.h"

#include <cmath>
#include <iostream>
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
  control::policy().filename() = "implosion_nx20.par";
  status = flecsi::start(control::execute);
  flecsi::finalize();
  return status;
}

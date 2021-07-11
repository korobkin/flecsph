#include "gtest/gtest.h"

#include <cmath>
#include <iostream>
#include "control.h"

#include <mpi.h>

namespace analysis {
enum e_conservation : size_t {
  MASS = 0,
  ENERGY = 1,
  MOMENTUM = 2,
  ANG_MOMENTUM = 3
};
}
using namespace analysis;
using namespace flecsi;

bool check_conservation(const std::vector<analysis::e_conservation> &);

int
main(int argc, char * argv[]) {
  auto status = flecsi::initialize(argc, argv);
  status = control::check_status(status);
  if(status != flecsi::run::status::success) {
    return status < flecsi::run::status::clean ? 0 : status;
  }
  control::policy().filename() = "KH_2d.par";
  flecsi::log::add_output_stream("clog", std::clog, true);
  assert(check_conservation({MASS, ENERGY, MOMENTUM}));
  status = flecsi::start(control::execute);
  flecsi::finalize();
  return status;
}

#include "gtest/gtest.h"
#include "flecsi/execution.hh"
#include <mpi.h>
#include <flecsi/flog.hh>
#include "control.h"

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
  control::policy().filename() = "sodtube_t1_n100.par";
  flecsi::log::add_output_stream("clog", std::clog, true);
  status = flecsi::start(control::execute);
  assert(check_conservation({MASS, ENERGY, MOMENTUM}));  
  flecsi::finalize();
  return status;
}
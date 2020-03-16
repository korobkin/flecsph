#include "gtest/gtest.h"

#include <cinchlog.h>
#include <iostream>
#include <cmath>

#include <mpi.h>

namespace flecsi{
namespace execution{
  void mpi_init_task(const char * parameter_file);
}
}

using namespace flecsi;
using namespace execution;

TEST(sedov, working) {
  MPI_Init(nullptr,nullptr);
  mpi_init_task("sedov_nx20.par");
  MPI_Finalize();
}

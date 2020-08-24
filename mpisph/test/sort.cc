#include "gtest/gtest.h"

#include <cmath>
#include <iostream>
#include <log.h>
#include <mpi.h>

#include "bodies_system.h"

using namespace ::testing;

// Rule for body equals == same position
inline bool
operator==(const body & b1, const body & b2) {
  return b1.coordinates() == b2.coordinates();
};

namespace flecsi {
namespace execution {
void
driver(int, char **) {}
} // namespace execution
} // namespace flecsi

TEST(tree_colorer, mpi_qsort) {
  MPI_Init(nullptr, nullptr);
  int rank;
  int size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  srand(time(NULL) * rank);
  log_set_output_rank(0);

  // Generating the particles randomly on each process
  int64_t nparticles = 100000;
  int64_t nparticlesperproc = nparticles / size;
  double maxbound = 1.0; // Particles positions between [0,1]
  // Adjust for last one
  if(rank == size - 1) {
    nparticlesperproc = (nparticles - nparticlesperproc * (size - 1));
  }
  log_one(info) << "Generating " << nparticles << std::endl;

  // Range to compute the keys
  std::array<point_t, 2> range;
  range[0] = point_t{};
  range[1] = point_t{maxbound, maxbound, maxbound};
  std::vector<body> bodies(nparticlesperproc);
  // Create the bodies and keys
  for(int64_t i = 0; i < nparticlesperproc; ++i) {
    // Random x, y and z
    bodies[i].set_coordinates(
      point_t{(double)rand() / (double)RAND_MAX * (maxbound),
        (double)rand() / (double)RAND_MAX * (maxbound),
        (double)rand() / (double)RAND_MAX * (maxbound)});

    // Compute the key
    bodies[i].set_key(key_type(range, bodies[i].coordinates()));
  }

  // Type used for sort 
  using sortType = std::pair<tree_topology_t::key_t,tree_topology_t::key_int_t>; 
  // Compare the sort type
  struct cmpType {
    bool operator()(const sortType& a, const sortType& b) const {
      if(a.first == b.first)
        return a.second < b.second; 
      return a.first < b.first; 
    }
  };
  struct extractType {
    sortType operator()(const body& a){
      return sortType(a.key(),a.id()); 
    }
  };

  auto bcomp = [](auto &left, auto &right) {
          if (left.key() < right.key()) {
            return true;
          }
          if (left.key() == right.key()) {
            return left.id() < right.id();
          }
          return false;
        };

  tree_colorer<body,sortType,extractType,cmpType> t; 
  t.hsort(bodies,nparticles, 
      bcomp); 

  // Check if the sort is valid: check if last particle of a rank 
  // is less than the first particle of next rank 
  assert(std::is_sorted(bodies.begin(), bodies.end(), bcomp));  

  using check_t = std::pair<key_type,key_type>; 

  check_t keys; 
  if(rank == 0){
    keys.first = key_type::min(); 
  }else{
    keys.first = bodies.front().key(); 
  }

  if(rank == size-1){
    keys.second = key_type::max(); 
  }else{
    keys.second = bodies.back().key(); 
  }

  std::vector<check_t> check(size); 

  MPI_Allgather(
    &keys, sizeof(check_t), MPI_BYTE,
    check.data(), sizeof(check_t), MPI_BYTE, MPI_COMM_WORLD); 

  MPI_Finalize();
}

#include "gtest/gtest.h"

#include <cmath>
#include <iostream>
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
  //int64_t nparticles = 1000000000;
  int64_t nparticlesperproc = 20000; //nparticles / size;
  int64_t nparticles = nparticlesperproc*size;
  double maxbound = 10000.0; // Particles positions between [0,1]
  // Adjust for last one
  if(rank == size - 1) {
    nparticlesperproc = (nparticles - nparticlesperproc * (size - 1));
  }
  log_one(info) << "Generating " << nparticles << " = " << nparticlesperproc <<" particles per process" << std::endl;

  MPI_Barrier(MPI_COMM_WORLD);
  auto start = omp_get_wtime();
  MPI_Barrier(MPI_COMM_WORLD);

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
  using sortType =
    std::pair<tree_topology_t::key_t, tree_topology_t::key_int_t>;
  // Compare the sort type
  struct cmpType {
    bool operator()(const sortType & a, const sortType & b) const {
      if(a.first == b.first)
        return a.second < b.second;
      return a.first < b.first;
    }
  };
  struct extractType {
    sortType operator()(const body & a) {
      return sortType(a.key(), a.id());
    }
  };
  struct cmpBody {
    bool operator()(const body & a, const body & b) const {
      if(a.key() == b.key())
        return a.id() < b.id();
      return a.key() < b.key();
    }
  };

  tree_colorer<sortType, body, extractType, cmpType, cmpBody> t;
  t.hsort(bodies, nparticles);

  MPI_Barrier(MPI_COMM_WORLD);
  if(rank ==  0)
  std::cout<<"time: "<<omp_get_wtime()-start<<std::endl;

  // Check if the sort is valid: check if last particle of a rank
  // is less than the first particle of next rank
  assert(std::is_sorted(bodies.begin(), bodies.end(), cmpBody{}));

  using check_t = std::pair<key_type, key_type>;

  check_t keys;

  keys.first = (rank == 0 ? key_type::min() : bodies.front().key());
  keys.second = (rank == size - 1 ? key_type::max() : bodies.back().key());

  std::vector<check_t> check(size);

  MPI_Allgather(&keys, sizeof(check_t), MPI_BYTE, check.data(), sizeof(check_t),
    MPI_BYTE, MPI_COMM_WORLD);

  if(rank == 0){
    // TODO: maybe std::is_sorted is better
    for(int i = 1 ; i < size ; ++i){
      if(!(check[i].first > check[i-1].second)){
        log_one(trace)<<i<<" ERROR: "<<check[i].first<<" !> "<<check[i-1].second<<std::endl;
      }
      assert(check[i].first > check[i-1].second);
    }
    // is this necessary?
    for(int i = 0 ; i < size-1 ; ++i){
      if(!(check[i].second < check[i+1].first)){
        log_one(trace)<<i<<" ERROR: "<<check[i].second<<" !< "<<check[i+1].first<<std::endl;
      }
      assert(check[i].second < check[i+1].first);

    }
  }


  MPI_Finalize();
}

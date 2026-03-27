#include "gtest/gtest.h"

#include <cmath>
#include <iostream>
#include <mpi.h>

#include "bodies_system.h"

// Rule for body equals == same position
inline bool
operator==(const body & b1, const body & b2) {
  return b1.coordinates() == b2.coordinates();
};



struct my_body{
  int64_t id() const {
    return 0; 
  }
  point_t coordinates() const {
    return coordinates_; 
  }
  key_type key() const {
    return key_; 
  }
  void set_coordinates(point_t p) { 
    coordinates_ = p; 
  } 
  void set_key(key_type k) {
    key_ = k; 
  }
  private: 
  point_t coordinates_;
  key_type key_;
}; 

int main(int argc, char* argv[]){
  MPI_Init(nullptr, nullptr);

  int64_t nparticles = 1000000;
  if(argc == 2 ){
    nparticles = atoll(argv[1]); 
  }else if(argc == 1){
    // Nothing 
  }else{
    std::cerr<<"ERROR ARGS"<<std::endl;
    MPI_Finalize(); 
    return 1; 
  }

  int rank;
  int size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  srand(time(NULL) * rank);
  log_set_output_rank(0);

  // Generating the particles randomly on each process
  int64_t nparticlesperproc = nparticles / size;
  double maxbound = 10000.0; // Particles positions between [0,1]
  // Adjust for last one
  if(rank == size - 1) {
    nparticlesperproc = (nparticles - nparticlesperproc * (size - 1));
  }
  log_one(info) << "Generating " << nparticles << " = " << nparticlesperproc <<" particles per process" << std::endl;
  log_one(info) << "Allocating: "<< 
    nparticlesperproc*sizeof(my_body)/1024./1024./1024. << " GB" 
    << " my_body: " << sizeof(my_body) <<" B"<< std::endl;

  // Range to compute the keys
  std::array<point_t, 2> range;
  range[0] = point_t{};
  range[1] = point_t{maxbound, maxbound, maxbound};
  std::vector<my_body> bodies(nparticlesperproc);
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
    sortType operator()(const my_body & a) {
      return sortType(a.key(), a.id());
    }
  };
  struct cmpBody {
    bool operator()(const my_body & a, const my_body & b) const {
      //if(a.key() == b.key())
      //  return a.id() < b.id();
      return a.key() < b.key();
    }
  };

  log_one(trace)<<"Sorting"<<std::endl;
  
  tree_colorer<sortType, my_body, extractType, cmpType, cmpBody> t;
  t.hsort(bodies, nparticles);

  log_one(trace)<<"Sort done"<<std::endl;;

  // Check if the sort is valid: check if last particle of a rank
  // is less than the first particle of next rank
  assert(std::is_sorted(bodies.begin(), bodies.end(), cmpBody{}));

  using check_t = std::pair<key_type, key_type>;

  check_t keys;
  if(rank == 0) {
    keys.first = key_type::min();
  }
  else {
    keys.first = bodies[0].key();
  }

  if(rank == size - 1) {
    keys.second = key_type::max();
  }
  else {
    keys.second = bodies.back().key();
  }

  std::vector<check_t> check(size);

  MPI_Allgather(&keys, sizeof(check_t), MPI_BYTE, check.data(), sizeof(check_t),
    MPI_BYTE, MPI_COMM_WORLD);

  if(rank == 0){
    for(int i = 1 ; i < size ; ++i){
      if(!(check[i].first > check[i-1].first)){
        log_one(trace)<<rank<<" ERROR: "<<check[rank].first<<" !< "<<check[rank-1].second<<std::endl;
      }
      assert(check[i].first > check[i-1].first); 
    }
    for(int i = 0 ; i < size-1 ; ++i){
      if(!(check[rank].second < check[rank+1].first)){
        log_one(trace)<<rank<<" ERROR: "<<check[rank].second<<" !< "<<check[rank+1].first<<std::endl;
      }
      assert(check[i].first > check[i-1].first); 

    }
  }
  MPI_Finalize();
}

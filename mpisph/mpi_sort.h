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
 * @file tree_colorer.h
 * @author Julien Loiseau
 * @date April 2017
 * @brief Function needed for MPI distribution of the bodies
 */

#pragma once

#include <fstream>
#include <iostream>
#include <numeric>
#include <omp.h>
#include <vector>

#include "default_physics.h"
#include "tree.h"
#include "utils.h"

#include "params.h" // For the variable smoothing length

using namespace mpi_utils;

// Output the data regarding the distribution for debug
#define OUTPUT_TREE_INFO 1

/**
 * @brief      All the function and buffers for the tree_colorer.
 *
 * @tparam     T     Type of the class
 * @tparam     D     Dimension for the problem
 * @todo fix it for type
 */
template<typename T> 
class tree_colorer {
private:

  using btype_t = T; 
  using splitter_t = std::pair<key_type, int64_t>; 
  using splitter_vector_t = std::vector<splitter_t>; 

  using histogram_t = std::vector<int>; 

  static constexpr auto splitter_t_comp = [](const splitter_t& a, const splitter_t& b)
  {
    if(a.first == b.first)
      return a.second < a.second; 
    return a.first < b.first; 
  }; 

  const int sample_factor_multiplier = 1; 
  const int sample_factor = 2; 
  const double epsilon_ = 0.000001; 
  const double max_rand = static_cast<double>(RAND_MAX); 

  MPI_Datatype MPI_T_SIZE_; 
  MPI_Datatype MPI_SPLITTER_SIZE_; 

  int size_; 
  int rank_;

  const int root_ = 0; 
  int k_ = 0; 
  int k = 0; 

  constexpr std::pair<int,int> target_range_(const int& N, const int& i, const int& p){
    return std::make_pair(
      (N*(i+1))/static_cast<double>(p)-(N*epsilon_)/(2.*p),
      (N*(i+1))/static_cast<double>(p)+(N*epsilon_)/(2.*p)); 
  }

  int max_sample_size_(){
    int size; 
    MPI_Comm_size(MPI_COMM_WORLD,&size); 
    return sample_factor*sample_factor_multiplier*size; 
  }

  // Can be changed it using openmp
  int sample_size_per_rank_(){
    return sample_factor*sample_factor_multiplier;
  }

  struct root_data {
    std::vector<splitter_vector_t> lower; 
    std::vector<splitter_vector_t> upper; 
  };

public:

  tree_colorer() {
    MPI_Comm_size(MPI_COMM_WORLD, &size_);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    // init mpi sizes 
    MPI_Type_contiguous(sizeof(btype_t), MPI_BYTE, &MPI_T_SIZE_);
    MPI_Type_commit(&MPI_T_SIZE_);

    MPI_Type_contiguous(sizeof(splitter_t), MPI_BYTE, &MPI_SPLITTER_SIZE_);
    MPI_Type_commit(&MPI_SPLITTER_SIZE_);

    k_ = log(log(size_)/epsilon_);
    srand(time(NULL));
  }

  ~tree_colorer() {}

  template<typename C> 
  void hsort(std::vector<btype_t> &rbodies, int totalnbodies, C&& comp ) {

    // 0. Sort local particles 
    std::sort(rbodies.begin(),rbodies.end(),comp); 

    int nprobe = sample_factor*sample_factor_multiplier; 
    const int total_nprobe = nprobe*size_; 

    const int nsplitters = size_ - 1;
    if(rank_ == root_){
      std::cout<<"k_ = "<<k_<<std::endl;
    }
    k = 0;  

    root_data rd; 
    if(rank_ == root_){
       // Allocate 
      rd.lower.resize(k_); 
      rd.upper.resize(k_); 
      for(int i = 0 ; i < k_ ; ++i){
        rd.lower[i].resize(nsplitters); 
        rd.upper[i].resize(nsplitters); 
      }
    }

    splitter_vector_t local_probe, global_probe;
    histogram_t hs, g_hs; 
    std::vector<int> lower_interval_rank(nsplitters); 
    std::vector<int> upper_interval_rank(nsplitters); 

    splitter_vector_t lower_interval(nsplitters,std::make_pair(key_type::min(),0)); 
    splitter_vector_t upper_interval(nsplitters,std::make_pair(key_type::max(),0));  

    do{

      sample_probe(totalnbodies,lower_interval,upper_interval,rbodies,nprobe,local_probe);  
      if(rank_ == root_){
        std::cout<<"sample_probe DONE"<<std::endl;
      }
      MPI_Barrier(MPI_COMM_WORLD); 
      gather_probe(local_probe,global_probe);
      if(rank_ == root_){
        std::cout<<"gather_probe DONE"<<std::endl;
      }
      MPI_Barrier(MPI_COMM_WORLD); 
      compute_histogram(global_probe,rbodies,hs); 
      if(rank_ == root_){
        std::cout<<"compute_histogram DONE"<<std::endl;
      }
      MPI_Barrier(MPI_COMM_WORLD); 
      reduce_histogram(hs,g_hs); 
      if(rank_ == root_){
        std::cout<<"reduce_histogram DONE"<<std::endl;
      }
      MPI_Barrier(MPI_COMM_WORLD); 
      // update root data 
      if(rank_ == root_){
        // Prefix sum the histogram 
        std::partial_sum(g_hs.begin(),g_hs.end(),g_hs.begin()); 
        assert(g_hs.back() == totalnbodies);
        // display histogram 
        std::cout<<std::endl<<"Histogram: "; 
        for(int i = 0 ; i < g_hs.size(); ++i){
          std::cout<<g_hs[i]<<"("<<i<<")"<<" - "; 
        }
        std::cout<<std::endl;
        std::cout<<std::endl<<"Objs: ";
        for(int i = 0 ; i < size_-1; ++i){
          auto rg = target_range_(totalnbodies, i, size_);
          std::cout<<(i+1)*totalnbodies/size_<<"("<<rg.first<<";"<<rg.second<<")"<<" - ";
        }
        std::cout<<std::endl;

        // Update L and U with values closest to the objectif  
        if(k > 0){
          std::copy(rd.lower[k-1].begin(),rd.lower[k-1].end(),rd.lower[k].begin());
          std::copy(rd.upper[k-1].begin(),rd.upper[k-1].end(),rd.upper[k].begin());
        }else{
          std::fill(rd.lower[k].begin(),rd.lower[k].end(),std::make_pair(key_type::min(),0));  
          std::fill(rd.upper[k].begin(),rd.upper[k].end(),std::make_pair(key_type::max(),0)); 
        }
        for(int i = 0 ; i < nsplitters ; ++i){ 
          // Range 
          int obj = totalnbodies*(i+1)/size_;
          for(int j = 0 ; j < g_hs.size(); ++j){
            // Lower boundary
            if(g_hs[j] <= obj){
              splitter_t splitter = (j==0) ? std::make_pair(key_type::min(),0L) : global_probe[j-1];
              if(splitter.first > rd.lower[k][i].first){
                rd.lower[k][i] = splitter; 
              }
            }
            // upper boundary
            if(g_hs[j] >= obj){
              splitter_t splitter = (j==g_hs.size() -1 )? std::make_pair(key_type::max(),0L) :global_probe[j];
              if(splitter.first < rd.upper[k][i].first)
                rd.upper[k][i] = splitter;
            }
          }
          assert(rd.lower[k][i] <= rd.upper[k][i]); 
        }
        //Choices 
        std::cout<<std::endl<<"Choices: "; 
        for(int i = 0 ; i < rd.lower[k].size(); ++i){
          std::cout<<rd.lower[k][i].first<<";"<<rd.upper[k][i].first<<" - "; 
        }
        // Change to 1 broadcast
        MPI_Bcast(rd.lower[k].data(), nsplitters, MPI_SPLITTER_SIZE_, root_, MPI_COMM_WORLD);
        MPI_Bcast(rd.upper[k].data(), nsplitters, MPI_SPLITTER_SIZE_, root_, MPI_COMM_WORLD);
        lower_interval = rd.lower[k]; 
        upper_interval = rd.upper[k]; 
      }else{
        // Change to 1 broadcast
        MPI_Bcast(lower_interval.data(), nsplitters, MPI_SPLITTER_SIZE_, root_, MPI_COMM_WORLD);
        MPI_Bcast(upper_interval.data(), nsplitters, MPI_SPLITTER_SIZE_, root_, MPI_COMM_WORLD);
      }

      MPI_Barrier(MPI_COMM_WORLD); 
      if(rank_ == root_){
        std::cout<<"Iteration: "<<k<<" DONE"<<std::endl;
      }
      ++k; 

    }while(k != k_); 

    // Take the middle of the interval and count the elements 
    splitter_vector_t final_splitters(size_-1);
    for(int i = 0 ; i < size_ -1; ++i){
      final_splitters[i] = std::make_pair(
        key_type(
          (static_cast<uint64_t>(lower_interval[i].first)+static_cast<uint64_t>(upper_interval[i].first))/2.),0L); 
    }
    // Reduction
    compute_histogram(final_splitters,rbodies,hs); 
    reduce_histogram(hs,g_hs); 
    if(rank_ == root_){
      std::cout<<"Result: ";
      std::cout<<" min("<<totalnbodies/size_-totalnbodies*epsilon_/2.*size_
        <<")-max("<<totalnbodies/size_+totalnbodies*epsilon_/2.*size_<<"): ";
      for(int i = 0 ; i < g_hs.size(); ++i){
        std::cout<<g_hs[i]<<" - ";
      }
      std::cout<<std::endl;
    }
    
    if(rank_ == root_)
      std::cout<<"DONE"<<std::endl;
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize(); 
    exit(0); 
  }


  void reduce_histogram(
    const histogram_t& hs, 
    histogram_t& g_hs)
  {
    if(rank_ == root_){
      g_hs.resize(hs.size());
      std::fill(g_hs.begin(),g_hs.end(),0);  
      assert(g_hs.size() == hs.size()); 
    }
    MPI_Reduce(
      hs.data(),g_hs.data(),hs.size(),MPI_INT,
      MPI_SUM,root_,MPI_COMM_WORLD);
  }

  void gather_probe(
    const splitter_vector_t& probe,
    splitter_vector_t& g_p)
  {
    // Send to all the number of probes
    int myprobes = probe.size(); 
    std::vector<int> nprobes(size_); 
    MPI_Allgather(&myprobes, 1, MPI_INT, nprobes.data(), 1,
     MPI_INT, MPI_COMM_WORLD); 
    std::vector<int> nprobes_displ(size_);

    std::partial_sum(nprobes.begin(), nprobes.end(), nprobes_displ.begin());
    nprobes_displ.insert(nprobes_displ.begin(),0);
    int total_nprobe = nprobes_displ.back(); 

    if(rank_ == root_)
      std::cout<<"Total probes: "<<total_nprobe<<std::endl;
    MPI_Barrier(MPI_COMM_WORLD); 

    g_p.resize(total_nprobe); 

    MPI_Allgatherv(probe.data(), probe.size(), MPI_SPLITTER_SIZE_,
                   g_p.data(), nprobes.data(), nprobes_displ.data(),
                   MPI_SPLITTER_SIZE_, MPI_COMM_WORLD);

    std::sort(g_p.begin(),g_p.end(),splitter_t_comp); 
  }

  // Generate sample in the interval 
  void sample_probe(
    const int64_t tnbodies, 
    const splitter_vector_t& lower_keys, 
    const splitter_vector_t& upper_keys, 
    const std::vector<btype_t>& bodies, 
    const int& nprobe, 
    splitter_vector_t& probe)
  {
    probe.clear(); 

    MPI_Barrier(MPI_COMM_WORLD);

    // Is it k+1 (k > 0) or k starts at 0? 
    double sampling_ratio = pow(2.*log(size_)/epsilon_,(k+1.)/static_cast<double>(k_));
    double proba = size_*sampling_ratio/static_cast<double>(tnbodies); 
    if(rank_ == root_){
      std::cout<<"Sampling_ratio: "<<sampling_ratio<<std::endl;
    }

    if(rank_ == root_){
      std::cout<<"Generating between: "<<std::endl; 
      for(int i = 0 ; i < lower_keys.size(); ++i){
        std::cout<<lower_keys[i].first<<" - "<<upper_keys[i].first<<std::endl;
      }
    }

    MPI_Barrier(MPI_COMM_WORLD);


    for(int i = 0 ; i < bodies.size(); ++i){
      bool find = false; 
      for(int j = 0; j < lower_keys.size(); ++j){
        if(bodies[i].key() > lower_keys[j].first && bodies[i].key() < upper_keys[j].first ){
          find = true; 
          break; 
        }
      }
      // Key is in range 
      if(find){
        double rnd = rand()/max_rand; 
        if(rnd < proba){
          probe.push_back(std::make_pair(bodies[i].key(),bodies[i].id()));
        }
      }
    }


    if(rank_ == root_){
      std::cout<<"Generated keys: "<<std::endl; 
      for(int i = 0 ; i < probe.size(); ++i){
        std::cout<<probe[i].first<<std::endl;
      }
    }

    MPI_Barrier(MPI_COMM_WORLD);
 
    std::cout<<rank_<<" nprobes = "<< probe.size() <<std::endl;
    MPI_Barrier(MPI_COMM_WORLD); 
  }

  // Compute the histogram and reduce histogram values 
  void compute_histogram(
    const splitter_vector_t& probe,
    const std::vector<btype_t>& bodies,
    histogram_t& hs)
  {
    hs.resize(probe.size()+1);
    std::fill(hs.begin(),hs.end(),0);  
    // 2. Compute histogram 
    int cur_probe = 0; 
    for(int i = 0; i < bodies.size(); ++i){
      if(bodies[i].key() < probe[cur_probe].first){
        ++hs[cur_probe]; 
      }else if(cur_probe == probe.size()){
        ++hs[cur_probe]; 
      }else{
        ++cur_probe;
        ++hs[cur_probe]; 
      }
    }
  }

}; // class tree_colorer

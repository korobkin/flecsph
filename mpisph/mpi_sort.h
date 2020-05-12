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

#ifndef _mpisph_tree_colorer_h_
#define _mpisph_tree_colorer_h_

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

  const size_t noct = 256 * 1024;   // Number of octets used for quicksort
  using btype_t = T; 
public:

  tree_colorer() {
    MPI_Comm_size(MPI_COMM_WORLD, &size_);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    // init mpi sizes 
    MPI_Type_contiguous(sizeof(btype_t), MPI_BYTE, &MPI_T_SIZE_);
    MPI_Type_commit(&MPI_T_SIZE_);

    MPI_Type_contiguous(sizeof(std::pair<key_type, int64_t>), MPI_BYTE, &MPI_PIVOT_SIZE_);
    MPI_Type_commit(&MPI_PIVOT_SIZE_);
  }

  ~tree_colorer() {}

  /*~---------------------------------------------------------------------------*
   * Function for sorting and distribution
   *~---------------------------------------------------------------------------*/

  /**
   * @brief      Sorting of the input particles or current particles using MPI.
   * This method is composed of several steps to implement the quick sort:
   * - Each process sorts its local particles
   * - Each process generate a subset of particles to fit the byte size limit
   * to send
   * - Each subset if send to the master (Here 0) who generates the pivot for
   * quick sort
   * - Pivot are send to each processes and they create buckets based on the
   * pivot
   * - Each bucket is send to the owner
   * - Each process sorts its local particles again
   * This is the first implementation, it can be long for the first sorting but
   * but then as the particles does not move very fast, the particles on the
   * edge are the only ones shared
   *
   * @param      rbodies       The rbodies, local bodies of this process.
   * @param[in]  totalnbodies  The totalnbodies on the overall simulation.
   */
  template<typename C> 
  void mpi_qsort(std::vector<btype_t> &rbodies, int totalnbodies,C&& comp ) {
    double s0 = omp_get_wtime(); 
    // Sort the keys
    // Use boost parallel sort
    std::sort(rbodies.begin(), rbodies.end(), comp); // sort

    // If one process, done
    if (size_ == 1) {
      log_one(trace) << "Local particles: " << totalnbodies << std::endl;
      return;
    } // if

    splitters_.clear();
    std::vector<int> scount(size_);
    generate_splitters_samples(splitters_, rbodies, totalnbodies);

    int cur_proc = 0;

    assert(splitters_.size() == size_ - 1 + 2);

    int64_t nbodies = rbodies.size();
    for (size_t i = 0L; i < nbodies; ++i) {
      if (rbodies[i].key() >= splitters_[cur_proc].first &&
          rbodies[i].key() < splitters_[cur_proc + 1].first) {
        scount[cur_proc]++;
      } else {
        i--;
        cur_proc++;
      }
    }

    // Check that we considered all the bodies
    assert(std::accumulate(scount.begin(), scount.end(), 0) == rbodies.size());

    std::vector<btype_t> recvbuffer;
    // Direct exchange using point to point
    mpi_alltoallv_p2p(scount, rbodies, recvbuffer);

    rbodies.clear();
    rbodies = recvbuffer;

    std::sort(rbodies.begin(), rbodies.end(), comp); // sort

    double s1 = omp_get_wtime(); 
#ifdef OUTPUT
    std::vector<int> totalprocbodies;
    totalprocbodies.resize(size_);
    int mybodies = rbodies.size();
    // Share the final array size of everybody
    MPI_Allgather(&mybodies, 1, MPI_INT, &totalprocbodies[0], 1, MPI_INT,
                  MPI_COMM_WORLD);
#ifdef OUTPUT_TREE_INFO
    std::ostringstream oss;
    oss << "Repartition (before): ";
    for (auto num : totalprocbodies)
      oss << num << ";";
    log_one(trace) << oss.str() <<"("<< s1-s0 <<"sec)"<< std::endl;
    oss.str("");
    oss.clear(); 
#endif
    //std::cout<<"BF = "<<rank_<<": "<<rbodies.size()<<": "<<rbodies.begin()->key()
    //  <<"("<<(rbodies.begin()+1)->key()<<") - ("<<(rbodies.end()-2)->key()<<")"<<rbodies.back().key()<<std::endl;
#endif // OUTPUT

MPI_Barrier(MPI_COMM_WORLD); 
double b0 = omp_get_wtime(); 
      balance_entities(rbodies, comp); 
double b1 = omp_get_wtime(); 

#ifdef OUTPUT
    mybodies = rbodies.size();
    // Share the final array size of everybody
    MPI_Allgather(&mybodies, 1, MPI_INT, &totalprocbodies[0], 1, MPI_INT,
                  MPI_COMM_WORLD);
#ifdef OUTPUT_TREE_INFO
    oss << "Repartition (After): ";
    for (auto num : totalprocbodies)
      oss << num << ";";
    log_one(trace) << oss.str() <<"("<< b1-b0 <<"sec)"<< std::endl;
#endif
    //std::cout<<"AF = "<<rank_<<": "<<rbodies.size()<<": "<<rbodies.begin()->key()
    //  <<"("<<(rbodies.begin()+1)->key()<<") - ("<<(rbodies.end()-2)->key()<<")"<<rbodies.back().key()<<std::endl;
#endif // OUTPUT
  }    // mpi_qsort

  template<typename C>
  void balance_entities(std::vector<btype_t>& data, C&&comp){

    MPI_Barrier(MPI_COMM_WORLD);

    // Gather what is on everyrank 
    std::vector<int> pcount(size_); 
    pcount[rank_] = data.size(); 
    MPI_Allgather(MPI_IN_PLACE, 1, MPI_INT, &pcount[0], 1, MPI_INT, MPI_COMM_WORLD);

    if(!rank_){
      int min_e = *std::min_element(pcount.begin(),pcount.end());
      int max_e = *std::max_element(pcount.begin(),pcount.end());
      log_one(trace)<<"Balance diff: ("<<min_e<<","<<max_e<<") = "<<max_e-min_e<<std::endl;
    }

    //if(rank == 0){
    //  std::cout<<"Count: -"; 
    //  for(int i = 0 ; i < size; ++i){
    //    std::cout<<pcount[i]<<"-";
    //  }
    //  std::cout<<std::endl;
    //}

    int64_t totalparticles = 0; 
    for(int i = 0 ; i < size_; ++i){
      totalparticles += pcount[i]; 
    }
    // Count the desire size
    std::vector<int> wdist(size_,0);
    int lparticles = totalparticles/size_; 
    int mparticles = totalparticles%size_; 

    for(int i = 0 ; i < size_; ++i){
      wdist[i] = lparticles; 
      if(i < mparticles){
        ++wdist[i]; 
      }
    }

    // Compute matrix of balance 
    std::vector<std::vector<int>> msend(size_); 
    std::vector<std::vector<int>> mrecv(size_);  

    for(int i = 0 ; i < size_ ; ++i){
      msend[i].resize(size_); 
      mrecv[i].resize(size_); 
    }

    for(int i = 0 ; i < size_-1 ; ++i){
      int need = pcount[i]-wdist[i];
      if(need < 0){
        need = -need; 
        // Need more
        assert(pcount[i+1] > need);
        mrecv[i][i+1] = need;
        msend[i+1][i] = need; 
        pcount[i+1] -= need;
        pcount[i] +=  need;
      }else if(need > 0){
        // Need less
        mrecv[i+1][i] = need;
        msend[i][i+1] = need;
        pcount[i] -= need;  
        pcount[i+1] += need; 
      }
    }

    //MPI_Barrier(MPI_COMM_WORLD);
    // Display matrix 
    //if(rank == 0){
    //  for(int i = 0 ; i < size; ++i){
    //    printf("%d(%05d) | ",i,wdist[i]); 
    //    for(int j = 0 ; j < size; ++j){
    //      printf("%05d(%05d)  ",mrecv[i][j],msend[i][j]); 
    //    }
    //    printf("\n"); 
    //  }
    //}
    //MPI_Barrier(MPI_COMM_WORLD);

    MPI_Request req; 

    // Do the send/receive to balance with Isend Irecv
    for(int i = 0 ; i < size_; ++i){
      if(msend[rank_][i] != 0){
        // Send lower or higher part of my particles 
        if(rank_ < i){
          MPI_Isend(
          &data[data.size()-msend[rank_][i]],
          msend[rank_][i],MPI_T_SIZE_,i,3,MPI_COMM_WORLD,
          &req); 
        }else{
          MPI_Isend(
          &data[0],
          msend[rank_][i],MPI_T_SIZE_,i,3,MPI_COMM_WORLD,
          &req); 
        }
      }
    }

    std::vector<btype_t> temp; 
    // Receive in a special buffer 
    for(int i = 0 ; i < size_ ; ++i){
      if(mrecv[rank_][i] != 0){
        temp.resize(temp.size()+mrecv[rank_][i]); 
        MPI_Recv(&temp[temp.size()-mrecv[rank_][i]],mrecv[rank_][i],MPI_T_SIZE_,i,3,MPI_COMM_WORLD,&Stat);
      }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // Clean local data and add temp 
    for(int i = 0 ; i < size_; ++i){
      if(msend[rank_][i] != 0){
        // Send lower or higher part of my particles 
        if(rank_ < i){
          data.erase(data.end()-msend[rank_][i],data.end());
        }else{
          data.erase(data.begin(),data.begin()+msend[rank_][i]);
        }
      }
    }

    data.insert(data.end(),temp.begin(),temp.end());
    std::sort(data.begin(),data.end(),comp);

    //localnparticles = data.size(); 
    //std::cout<<rank<<"SDS_SORT done: "<<data.size()<<": "<<data.begin()->key()
    //  <<"("<<(data.begin()+1)->key()<<") - ("<<(data.end()-2)->key()<<")"<<data.back().key()<<std::endl;
    
    MPI_Barrier(MPI_COMM_WORLD); 
  }

  /**
   * @brief      Use in mpi_qsort to generate the splitters to sort the
   * particles in the quick sort algorithm In this function we take some
   * samplers of the total particles and the root determines the splitters This
   * version is based on the sample splitter algorithm but we generate more
   * samples on each process
   *
   * @param      splitters  The splitters used in the qsort in mpi_qsort
   * @param[in]  rbodies  The local bodies of the process
   */
  void generate_splitters_samples(
      std::vector<std::pair<key_type, int64_t>> &splitters,
      std::vector<btype_t> &rbodies, const int64_t totalnbodies) {

    // Create a vector for the samplers
    std::vector<std::pair<key_type, int64_t>> keys_sample;
    // Number of elements for sampling
    // In this implementation we share up to 256KB to
    // the master.
    size_t maxnsamples = noct / sizeof(std::pair<key_type, int64_t>);
    int64_t nvalues = rbodies.size();
    size_t nsample = maxnsamples * ((double)nvalues / (double)totalnbodies);

    if (nvalues < (int64_t)nsample) {
      nsample = nvalues;
    }

    for (size_t i = 0; i < nsample; ++i) {
      int64_t position = (nvalues / (nsample + 1.)) * (i + 1.);
      keys_sample.push_back(
          std::make_pair(rbodies[position].key(), rbodies[position].id()));
    } // for
    assert(keys_sample.size() == (size_t)nsample);

    std::vector<std::pair<key_type, int64_t>> master_keys;
    std::vector<int> master_recvcounts;
    std::vector<int> master_offsets;
    int master_nkeys = 0;

    if (rank_ == 0) {
      master_recvcounts.resize(size_);
    } // if

    // Echange the number of samples
    MPI_Gather(&nsample, 1, MPI_INT, &master_recvcounts[0], 1, MPI_INT, 0,
               MPI_COMM_WORLD);

    // Master
    // Sort the received keys and create the pivots
    if (rank_ == 0) {
      master_offsets.resize(size_);
      master_nkeys = std::accumulate(master_recvcounts.begin(),
                                     master_recvcounts.end(), 0);
      if (totalnbodies < master_nkeys) {
        master_nkeys = totalnbodies;
      }
      // Number to receiv from each process
      for (int i = 0; i < size_; ++i) {
        master_recvcounts[i] *= sizeof(std::pair<key_type, int64_t>);
      } // for
      std::partial_sum(master_recvcounts.begin(), master_recvcounts.end(),
                       &master_offsets[0]);
      master_offsets.insert(master_offsets.begin(), 0);
      master_keys.resize(master_nkeys);
    } // if

    MPI_Gatherv(&keys_sample[0], nsample * sizeof(std::pair<key_type, int64_t>),
                MPI_BYTE, &master_keys[0], &master_recvcounts[0],
                &master_offsets[0], MPI_BYTE, 0, MPI_COMM_WORLD);

    // Generate the splitters, add zero and max keys
    splitters.resize(size_ - 1 + 2);
    if (rank_ == 0) {
      std::sort(master_keys.begin(), master_keys.end(),
                [](auto &left, auto &right) {
                  if (left.first < right.first) {
                    return true;
                  }
                  if (left.first == right.first) {
                    return left.second < right.second;
                  }
                  return false;
                });

      splitters[0].first = key_type::min();
      splitters[0].second = 0L;
      splitters[size_].first = key_type::max();
      splitters[size_].second = LONG_MAX;

      for (int i = 0; i < size_ - 1; ++i) {
        int64_t position = (master_nkeys / size_) * (i + 1);
        splitters[i + 1] = master_keys[position];
        assert(splitters[i + 1].first > splitters[0].first &&
               splitters[i + 1].first < splitters[size_].first);
      } // for
    } // if

    // Bradcast the splitters
    MPI_Bcast(&splitters[0],
              (size_ - 1 + 2) * sizeof(std::pair<key_type, int64_t>), MPI_BYTE,
              0, MPI_COMM_WORLD);

    if(!rank_){
      std::ostringstream oss;
      oss << "Pivots: ";
      for(int i = 0 ; i < splitters.size(); ++i){
        oss<<splitters[i].first<<": ";
      } 
      log_one(trace)<<oss.str()<<std::endl;
    }
  }

private:
  // Key track of the splitter to know first and last key
  std::vector<std::pair<key_type, int64_t>> splitters_;

  MPI_Datatype MPI_T_SIZE_; 
  MPI_Datatype MPI_PIVOT_SIZE_; 

  int size_, rank_;

}; // class tree_colorer

#endif // _mpisph_tree_colorer_h_

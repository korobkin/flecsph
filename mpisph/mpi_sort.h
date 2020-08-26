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

/**
 * @brief      Histogram distributed sort.
 * Implementation based on: https://arxiv.org/pdf/1803.01237.pdf
 *
 * @tparam Key          Type of keys for sorting
 * @tparam Type         Type of entities sorted
 * @tparam Extract      Function to extract the key from the type
 * @tparam CompareKey   Function to compare two Key types
 * @tparam COmpareType  Function two compare two Types
 */
template<
  typename Key,
  typename Type,
  typename Extract,
  typename Compare = std::less<Key>,
  typename CompareType = std::less<Type>>
class tree_colorer
{

public:
  using btype_t = Type;
  using splitter_t = Key;
  using splitter_vector_t = std::vector<splitter_t>;
  using histogram_t = std::vector<int>;
  using compare_t = Compare;
  using compare_type_t = CompareType; 
  using extract_t = Extract;

  tree_colorer() {
    MPI_Comm_size(MPI_COMM_WORLD, &size_);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    // init mpi sizes
    MPI_Type_contiguous(sizeof(btype_t), MPI_BYTE, &MPI_T_SIZE_);
    MPI_Type_commit(&MPI_T_SIZE_);

    MPI_Type_contiguous(sizeof(splitter_t), MPI_BYTE, &MPI_SPLITTER_SIZE_);
    MPI_Type_commit(&MPI_SPLITTER_SIZE_);

    nsplitters_ = size_ - 1;
  }

  ~tree_colorer() {}

  void hsort(std::vector<btype_t> & rbodies,
    int totalnbodies,
    const double epsilon = 0.05) {

    srand(time(NULL) * rank_);

    epsilon_ = epsilon;
    k_ = log(log(size_) / epsilon_);

    std::sort(rbodies.begin(), rbodies.end(), compare_type);

    if(size_ == 1) {
      return;
    }

    log_one(trace) << "k_ = " << k_ << " nsplitters = " << nsplitters_
                << std::endl;
    k = 0;

    // Root only data
    std::vector<splitter_vector_t> lower;
    std::vector<splitter_vector_t> upper;
    if(rank_ == root_) {
      // Allocate
      lower.resize(k_);
      upper.resize(k_);
      for(int i = 0; i < k_; ++i) {
        lower[i].resize(nsplitters_);
        upper[i].resize(nsplitters_);
      }
    }

    splitter_vector_t probes;
    histogram_t hs;
    std::vector<int> lower_interval_rank(nsplitters_);
    std::vector<int> upper_interval_rank(nsplitters_);

    splitter_vector_t lower_interval(
      nsplitters_, splitter_t(key_type::min(), 0));
    splitter_vector_t upper_interval(
      nsplitters_, splitter_t(key_type::max(), 0));

    do {

      sample_allgather_probe_(
        totalnbodies, lower_interval, upper_interval, rbodies, probes);
      compute_reduce_histogram_(probes, rbodies, hs);
      // update root data
      if(rank_ == root_) {
        // Prefix sum the histogram
        std::partial_sum(hs.begin(), hs.end(), hs.begin());
        assert(hs.back() == totalnbodies);
        // display histogram
        //std::cout << std::endl << "Histogram: ";
        //for(int i = 0; i < hs.size(); ++i) {
        //  std::cout << hs[i] << "(" << i << ")"
        //            << " - ";
        //}
        //std::cout << std::endl;
        //std::cout << std::endl << "Objs: ";
        //for(int i = 0; i < size_ - 1; ++i) {
        //  auto rg = target_range_(totalnbodies, i, size_);
        //  std::cout << (i + 1) * totalnbodies / size_ << "(" << rg.first << ";"
        //            << rg.second << ")"
        //            << " - ";
        //}
        //std::cout << std::endl;

        // Update L and U with values closest to the objectif
        if(k > 0) {
          std::copy(lower[k - 1].begin(), lower[k - 1].end(), lower[k].begin());
          std::copy(upper[k - 1].begin(), upper[k - 1].end(), upper[k].begin());
        }
        else {
          std::fill(
            lower[k].begin(), lower[k].end(), splitter_t(key_type::min(), 0));
          std::fill(
            upper[k].begin(), upper[k].end(), splitter_t(key_type::max(), 0));
        }
        for(int i = 0; i < nsplitters_; ++i) {
          if(upper[k][i].first == lower[k][i].first)
            continue; 
          // Range
          int obj = totalnbodies * (i + 1) / size_;
          for(int j = 0; j < hs.size(); ++j) {
            // Lower boundary
            if(hs[j] <= obj) {
              splitter_t splitter =
                (j == 0) ? splitter_t(key_type::min(), 0) : probes[j - 1];
              if(splitter.first > lower[k][i].first) {
                lower[k][i] = splitter;
              }
            }
            // upper boundary
            if(hs[j] >= obj) {
              splitter_t splitter = (j == hs.size() - 1)
                                      ? splitter_t(key_type::max(), 0)
                                      : probes[j];
              if(splitter.first < upper[k][i].first)
                upper[k][i] = splitter;
            }
          }
        }
        // Choices
        //std::cout << std::endl << "Choices: ";
        //for(int i = 0; i < lower[k].size(); ++i) {
        //  std::cout << lower[k][i].first << ";" << upper[k][i].first << " - ";
        //}
        for(int i = 0; i < lower[k].size(); ++i) {
          assert(lower[k][i] <= upper[k][i]);
        }
        // Change to 1 broadcast
        MPI_Bcast(lower[k].data(), nsplitters_, MPI_SPLITTER_SIZE_, root_,
          MPI_COMM_WORLD);
        MPI_Bcast(upper[k].data(), nsplitters_, MPI_SPLITTER_SIZE_, root_,
          MPI_COMM_WORLD);
        lower_interval = lower[k];
        upper_interval = upper[k];
      }
      else {
        // Change to 1 broadcast
        MPI_Bcast(lower_interval.data(), nsplitters_, MPI_SPLITTER_SIZE_, root_,
          MPI_COMM_WORLD);
        MPI_Bcast(upper_interval.data(), nsplitters_, MPI_SPLITTER_SIZE_, root_,
          MPI_COMM_WORLD);
      }

      MPI_Barrier(MPI_COMM_WORLD);
      //if(rank_ == root_) {
      //  std::cout << "Iteration: " << k << " DONE" << std::endl;
      //}
      ++k;

    } while(k < k_);

    // Latest histogram
    //if(rank_ == root_) {
    //  std::cout << std::endl << "Lastest Histogram: ";
    //  for(int i = 0; i < hs.size(); ++i) {
    //    std::cout << hs[i] << "(" << i << ")"
    //              << " - ";
    //  }
    //  std::cout << std::endl;
    //}

    // Take the middle of the interval and count the elements
    splitter_vector_t final_splitters(size_ - 1);
    final_splitters = upper_interval;

    //if(rank_ == root_) {
    //  std::cout << std::endl << "Final splitters: ";
    //  for(int i = 0; i < size_ - 1; ++i) {
    //    std::cout << final_splitters[i].first << " - ";
    //  }
    //  std::cout << std::endl;
    //}

    // Reduction
    compute_reduce_histogram_(final_splitters, rbodies, hs);
    if(rank_ == root_) {
      auto rg = target_range_(totalnbodies, 0, size_);
      std::ostringstream oss; 

      oss << "Splitters: ";
      oss << " [" << rg.first << ";" << rg.second << "]: ";
      for(int i = 0; i < hs.size(); ++i) {
        oss << hs[i];
        if(!(hs[i] >= rg.first && hs[i] <= rg.second)) {
          oss << ":F";
        }
        oss << " - ";
      }
      log_one(trace) << oss.str() << std::endl;
      //for(int i = 0 ; i < hs.size(); ++i){
      //  assert(hs[i] >= rg.first && hs[i] <= rg.second);
      //}
    }

    // Use the splitters to distribute data
    exchange_entities_(final_splitters, rbodies);

  }

private:

  void exchange_entities_(splitter_vector_t & splitters,
    std::vector<btype_t> & bodies) 
  {
    // Generate local buckets
    std::vector<int> offsets(size_);
    int cur_splitter = 0;
    for(int i = 0; i < bodies.size(); ++i) {
      if(compare_key(extract(bodies[i]),splitters[cur_splitter])) {
        ++(offsets[cur_splitter]);
      }
      else if(cur_splitter == size_ - 1) {
        ++(offsets[cur_splitter]);
      } else {
        --i; 
        ++cur_splitter; 
      }
    }

    std::vector<int> recvcount(size_), recvoffsets(size_), sendoffsets(size_);
    // Exchange the send count
    MPI_Alltoall(
      &offsets[0], 1, MPI_INT, &recvcount[0], 1, MPI_INT, MPI_COMM_WORLD);

    std::partial_sum(recvcount.begin(), recvcount.end(), recvoffsets.begin());
    recvoffsets.insert(recvoffsets.begin(), 0);
    std::partial_sum(offsets.begin(), offsets.end(), sendoffsets.begin());
    sendoffsets.insert(sendoffsets.begin(), 0);
    // Set the recvbuffer to the right size
    std::vector<btype_t> recvbuffer;
    recvbuffer.resize(recvoffsets.back());

    std::vector<MPI_Status> status(size_);
    std::vector<MPI_Request> request(size_);
    for(int i = 0; i < size_; ++i) {
      if(offsets[i] != 0) {
        auto * start = bodies.data();
        MPI_Isend(start + sendoffsets[i], offsets[i], MPI_T_SIZE_, i, 0,
          MPI_COMM_WORLD, &request[i]);
      }
    }
    for(int i = 0; i < size_; ++i) {
      if(recvcount[i] != 0) {
        auto * start = recvbuffer.data();
        MPI_Recv(start + recvoffsets[i], recvcount[i], MPI_T_SIZE_, i, MPI_ANY_TAG,
          MPI_COMM_WORLD, &status[i]);
      }
      if(offsets[i] != 0) {
        MPI_Wait(&request[i], &status[i]);
      }
    }

    bodies = recvbuffer;

    // Sort end buffer
    std::sort(bodies.begin(), bodies.end(), compare_type);
  }

  // Generate sample in the interval
  void sample_allgather_probe_(const int64_t tnbodies,
    const splitter_vector_t & lower_keys,
    const splitter_vector_t & upper_keys,
    const std::vector<btype_t> & bodies,
    splitter_vector_t & probes) {
    splitter_vector_t local_probes;

    // Is it k+1 (k > 0) or k starts at 0?
    const double sampling_ratio =
      pow(2. * log(size_) / epsilon_, (k + 1.) / static_cast<double>(k_));
    const double proba = size_ * sampling_ratio / static_cast<double>(tnbodies);
    
    for(int i = 0; i < bodies.size(); ++i) {
      bool find = false;
      for(int j = 0; j < lower_keys.size(); ++j) {
        if(lower_keys[j] != upper_keys[j]){
          if(compare_key(lower_keys[j],extract(bodies[i]))  &&
            compare_key(extract(bodies[i]),upper_keys[j])) {
            find = true;
            break;
          }
        }
      }
      // Key is in range
      if(find) {
        double rnd = rand() / max_rand;
        if(rnd < proba) {
          local_probes.push_back(splitter_t(extract(bodies[i])));
        }
      }
    }

    // If first iteration, force at least one probe per rank
    if( k == 0 && local_probes.size() == 0) {
      int middle = bodies.size() / 2;
      local_probes.push_back(
        splitter_t(extract(bodies[middle])));
    }

    // Send to all the number of probes
    int myprobes = local_probes.size();
    std::vector<int> nprobes(size_);
    MPI_Allgather(
      &myprobes, 1, MPI_INT, nprobes.data(), 1, MPI_INT, MPI_COMM_WORLD);
    std::vector<int> nprobes_displ(size_);

    std::partial_sum(nprobes.begin(), nprobes.end(), nprobes_displ.begin());
    nprobes_displ.insert(nprobes_displ.begin(), 0);
    int total_nprobe = nprobes_displ.back();

    probes.resize(total_nprobe);

    MPI_Allgatherv(local_probes.data(), myprobes, MPI_SPLITTER_SIZE_,
      probes.data(), nprobes.data(), nprobes_displ.data(), MPI_SPLITTER_SIZE_,
      MPI_COMM_WORLD);

    std::sort(probes.begin(), probes.end(), compare_key);

  }

  // Compute the histogram and reduce histogram values
  void compute_reduce_histogram_(const splitter_vector_t & probe,
    const std::vector<btype_t> & bodies,
    histogram_t & hs) {
    hs.resize(probe.size() + 1);
    std::fill(hs.begin(), hs.end(), 0);
    // 2. Compute histogram
    int cur_probe = 0;
    for(int i = 0; i < bodies.size(); ++i) {
      if(compare_key(extract(bodies[i]),probe[cur_probe])) {
        ++(hs[cur_probe]);
      }
      else if(cur_probe == probe.size()) {
        ++(hs[cur_probe]);
      }
      else {
        --i; 
        ++cur_probe;
      }
    }
    if(rank_ == root_) {
      MPI_Reduce(MPI_IN_PLACE, hs.data(), hs.size(), MPI_INT, MPI_SUM, root_,
        MPI_COMM_WORLD);
    }
    else {
      MPI_Reduce(
        hs.data(), nullptr, hs.size(), MPI_INT, MPI_SUM, root_, MPI_COMM_WORLD);
    }
  }

  compare_t compare_key;
  compare_type_t compare_type; 
  extract_t extract;

  // Percent of keys around splitters
  double epsilon_ = 0.05;
  const double max_rand = static_cast<double>(RAND_MAX);
  // Rank being root for reduce
  const int root_ = 0;

  MPI_Datatype MPI_T_SIZE_;
  MPI_Datatype MPI_SPLITTER_SIZE_;

  // MPI data
  int size_, rank_;
  // k_ number of iterations, k current iteration
  int k_ = 0, k = 0;
  // Total number of splitters to generate
  int nsplitters_ = 0;

  constexpr std::pair<int, int>
  target_range_(const int & N, const int & i, const int & p) {
    return std::make_pair(
      (N * (i + 1)) / static_cast<double>(p) - (N * epsilon_) / (2. * p),
      (N * (i + 1)) / static_cast<double>(p) + (N * epsilon_) / (2. * p));
  }

}; // class tree_colorer

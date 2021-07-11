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
#include <vector>
#include <random>

#include "omp.h"
#include "mpi.h"

constexpr double FLECSPH_HSORT_DEFAULT_EPS = 0.05;
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
  using histogram_t = std::vector<int64_t>;
  using compare_t = Compare;
  using compare_type_t = CompareType;
  using extract_t = Extract;

  using interval_t = std::pair<splitter_t, splitter_t>;
  using interval_vector_t = std::vector<interval_t>;

  tree_colorer() :
    rng_gen_(rng_dev_()),
    rng_dist_(0.0, 1.0) {
    MPI_Comm_size(MPI_COMM_WORLD, &size_);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    // init mpi sizes
    MPI_Type_contiguous(sizeof(btype_t), MPI_BYTE, &MPI_T_SIZE_);
    MPI_Type_commit(&MPI_T_SIZE_);

    MPI_Type_contiguous(sizeof(splitter_t), MPI_BYTE, &MPI_SPLITTER_SIZE_);
    MPI_Type_commit(&MPI_SPLITTER_SIZE_);

    MPI_Type_contiguous(sizeof(interval_t), MPI_BYTE, &MPI_INTERVAL_SIZE_);
    MPI_Type_commit(&MPI_INTERVAL_SIZE_);

    // number of particle splitters
    nsplitters_ = size_ - 1;
    // the "interval" containing the splitter
    intervals_.resize(nsplitters_);
    // set this as default
    setEpsilon(FLECSPH_HSORT_DEFAULT_EPS);
  }

  ~tree_colorer() {
    // free the data-type allocations
    // need to reorder sort test for this
    /*MPI_Type_free(&MPI_T_SIZE_);
    MPI_Type_free(&MPI_SPLITTER_SIZE_);
    MPI_Type_free(&MPI_INTERVAL_SIZE_);*/
  }

  // get/set for the sort epsilon
  // this lets us evaluate some vars once
  inline void setEpsilon(const double eps) {
    epsilon_ = eps;
    // theorem 4.8;
    // NOTE: C rounds -> zero (truncates), so add 0.5 to do normal rounding
    nrounds_ = static_cast<int>(std::log(std::log(size_)/epsilon_) + 0.5);
  }
  constexpr inline auto& getEpsilon() const { return epsilon_; }

  void hsort(std::vector<btype_t> & rbodies,
    int64_t totalnbodies)
  {

    std::sort(rbodies.begin(), rbodies.end(), compare_type);

    if(size_ == 1) { return; }

    log_one(trace) << "nrounds = " << nrounds_ << " nsplitters = " << nsplitters_
                   << std::endl;

    splitter_vector_t probes;
    histogram_t hs;

    std::fill(std::begin(intervals_), std::end(intervals_), std::make_pair(splitter_t(key_type::min(), 0), splitter_t(key_type::max(),0)));

    for(int hitr=0; hitr < nrounds_; ++hitr) {
      log_one(trace) << "hitr: " << hitr << std::endl;

      sample_allgather_probe_(
        totalnbodies, rbodies, hitr, probes);
      compute_reduce_histogram_(probes, rbodies, hs);
      // update root data
      if(rank_ == root_) {
        // Prefix sum the histogram
        std::partial_sum(hs.begin(), hs.end(), hs.begin());
        assert(hs.back() == totalnbodies);
        // Update L and U with values closest to the objectif
        for(int i = 0; i < nsplitters_; ++i) {
          auto& [lo_s, hi_s] = intervals_[i];
          if(lo_s.first == hi_s.first)
            continue;
          // Range
          int64_t obj = totalnbodies * (i + 1) / size_;
          for(int j = 0; j < hs.size(); ++j) {
            // Lower boundary
            if(hs[j] <= obj) {
              splitter_t splitter =
                (j == 0) ? splitter_t(key_type::min(), 0) : probes[j - 1];
              if(splitter.first > lo_s.first) {
                lo_s = splitter;
              }
            }
            // upper boundary
            if(hs[j] >= obj) {
              splitter_t splitter = (j == hs.size() - 1)
                                      ? splitter_t(key_type::max(), 0)
                                      : probes[j];
              if(splitter.first < hi_s.first)
                hi_s = splitter;
            }
          }
        }
        for(const auto& [lo_s, hi_s]: intervals_)
          assert(lo_s <= hi_s);
      }
      MPI_Bcast(intervals_.data(), nsplitters_, MPI_INTERVAL_SIZE_, root_, MPI_COMM_WORLD);

      MPI_Barrier(MPI_COMM_WORLD);
    }
    // Take the middle of the interval and count the elements
    splitter_vector_t final_splitters;
    for(const auto& [lo_s, hi_s]: intervals_)
    {
      final_splitters.push_back(hi_s);
    }

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
    std::vector<int64_t> offsets(size_);
    int cur_splitter = 0;
    for(int64_t i = 0; i < bodies.size(); ++i) {
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

    std::vector<int64_t> recvcount(size_), recvoffsets(size_), sendoffsets(size_);
    // Exchange the send count
    MPI_Alltoall(
      &offsets[0], 1, MPI_INT64_T, &recvcount[0], 1, MPI_INT64_T, MPI_COMM_WORLD);

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
    const std::vector<btype_t> & bodies,
    const int round,
    splitter_vector_t& probes) {

    splitter_vector_t sample_space, local_probes;

    // Is it k+1 (k > 0) or k starts at 0?
    const double sampling_ratio =
      pow(2. * std::log(size_) / epsilon_, (static_cast<double>(round) + 1.) / static_cast<double>(nrounds_));
    const double proba = size_ * sampling_ratio / static_cast<double>(tnbodies);

    // first, generate the sample space with keys within intervals
    for(const auto& bod : bodies) {
      for(auto&& [lo_k, up_k] : intervals_)
      {
        if(lo_k != up_k)
        {
          auto bkey = extract(bod);
          if(compare_key(lo_k, bkey) && compare_key(bkey, up_k))
          {
            sample_space.push_back(splitter_t(bkey));
            break;
          }
        }
      }
    }
    // second, select probes by sampling the sample space
    std::copy_if(std::begin(sample_space), std::end(sample_space), std::back_inserter(local_probes), [&](auto&&){ return (rng_dist_(rng_gen_) < proba); });

    // If first iteration, force at least one probe per rank
    if( round == 0 && local_probes.size() == 0) {
      int64_t middle = bodies.size() / 2;
      local_probes.push_back(
        splitter_t(extract(bodies[middle])));
    }

    // Send to all the number of probes
    int n_local_probes = local_probes.size();
    std::vector<int> nprobes(size_);
    MPI_Allgather(
      &n_local_probes, 1, MPI_INT, nprobes.data(), 1, MPI_INT, MPI_COMM_WORLD);
    std::vector<int> nprobes_displ(size_);

    std::partial_sum(nprobes.begin(), nprobes.end(), nprobes_displ.begin());
    nprobes_displ.insert(nprobes_displ.begin(), 0);
    int total_nprobe = nprobes_displ.back();

    probes.resize(total_nprobe);

    MPI_Allgatherv(local_probes.data(), n_local_probes, MPI_SPLITTER_SIZE_,
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
    for(int64_t i = 0; i < bodies.size(); ++i) {
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
      MPI_Reduce(MPI_IN_PLACE, hs.data(), hs.size(), MPI_INT64_T, MPI_SUM, root_,
        MPI_COMM_WORLD);
    }
    else {
      MPI_Reduce(
        hs.data(), nullptr, hs.size(), MPI_INT64_T, MPI_SUM, root_, MPI_COMM_WORLD);
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
  MPI_Datatype MPI_INTERVAL_SIZE_;

  // MPI data
  int size_, rank_;
  // k_ number of iterations, k current iteration
  int nrounds_ = 0;
  // Total number of splitters to generate
  int nsplitters_ = 0;

  interval_vector_t intervals_;

  // for RNG
  std::random_device rng_dev_;
  std::mt19937 rng_gen_;
  std::uniform_real_distribution<> rng_dist_;

  constexpr std::pair<int, int>
  target_range_(const int & N, const int & i, const int & p) {
    return std::make_pair(
      (N * (i + 1)) / static_cast<double>(p) - (N * epsilon_) / (2. * p),
      (N * (i + 1)) / static_cast<double>(p) + (N * epsilon_) / (2. * p));
  }

}; // class tree_colorer

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

#ifdef BOOST
#include <boost/sort/sort.hpp>
#endif 

#include "default_physics.h"
#include "tree.h"
#include "utils.h"

#include "params.h" // For the variable smoothing length

using namespace mpi_utils;

//#define BOOST_PARALLEL 1
// Output the data regarding the distribution for debug
#define OUTPUT_TREE_INFO 1

/**
 * Structrue for branch distribution
 */
struct mpi_branch_t {
  point_t coordinates;
  double mass;
  point_t min;
  point_t max;
  key_type key;
  int owner;
  size_t sub_entities;
  bool leaf;
};

/**
 * @brief      All the function and buffers for the tree_colorer.
 *
 * @tparam     T     Type of the class
 * @tparam     D     Dimension for the problem
 * @todo fix it for type
 */
template <typename T, size_t D> class tree_colorer {
private:
  const size_t noct = 1024 * 1024;   // Number of octets used for quicksort

public:
  static const size_t dimension = D;
  using point_t = flecsi::point_u<T, dimension>;

  tree_colorer() {}

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
  void mpi_qsort(std::vector<body> &rbodies, int totalnbodies) {
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    double timer = omp_get_wtime(); 

    std::sort(rbodies.begin(), rbodies.end(), 
      [](auto &left, auto &right) {
          if (left.key() < right.key()) {
            return true;
          }
          if (left.key() == right.key()) {
            return left.id() < right.id();
          }
          return false;
      }); // sort

    if (size == 1) {
      clog_one(trace) << "Local particles: " << totalnbodies << std::endl;
      return;
    } // if

    splitters_.clear();
    std::vector<int> scount(size);
    generate_splitters_samples(splitters_, rbodies, totalnbodies);

    int cur_proc = 0;

    assert(splitters_.size() == size - 1 + 2);

    int64_t nbodies = rbodies.size();
    for (size_t i = 0L; i < nbodies; ++i) {
      if (rbodies[i].key() >= splitters_[cur_proc].first &&
          rbodies[i].key() < splitters_[cur_proc + 1].first) {
        scount[cur_proc]++;
      } else {
        i--;
        cur_proc++;
      } // if
    } // for

    // Check that we considered all the bodies
    assert(std::accumulate(scount.begin(), scount.end(), 0) == rbodies.size());

    std::vector<body> recvbuffer;
    // Direct exchange using point to point
    mpi_alltoallv_p2p(scount, rbodies, recvbuffer);

    rbodies.clear();
    rbodies = recvbuffer;

    std::sort(rbodies.begin(), rbodies.end(), 
      [](auto &left, auto &right) {
          if (left.key() < right.key()) {
            return true;
          }
          if (left.key() == right.key()) {
            return left.id() < right.id();
          }
          return false;
      }); // sort

#ifdef OUTPUT_TREE_INFO
    std::vector<int> totalprocbodies;
    totalprocbodies.resize(size);
    int mybodies = rbodies.size();
    // Share the final array size of everybody
    MPI_Allgather(&mybodies, 1, MPI_INT, &totalprocbodies[0], 1, MPI_INT,
                  MPI_COMM_WORLD);
    int min = *std::min_element(totalprocbodies.begin(),totalprocbodies.end()); 
    int max = *std::max_element(totalprocbodies.begin(),totalprocbodies.end());
    clog_one(trace) <<std::fixed<<std::setprecision(2)<< "Repartition: min("<<
      min<<") max("<<max<<") diff="<<max-min<<" "<<
      omp_get_wtime()-timer<<"s"<<std::endl;
#endif // OUTPUT_TREE_INFO
  } // mpi_qsort

  /*~---------------------------------------------------------------------------*
   * Utils functions
   *~---------------------------------------------------------------------------*/

  /**
   * @brief      Compute the global range of all the particle system
   *
   * @param      bodies           The bodies, local of this process
   * @param      range            The range computed in that function
   */
  void mpi_compute_range(const std::vector<body> &bodies,
                         std::array<point_t, 2> &range) {
    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Compute the local range
    range_t lrange;

    lrange[1] = bodies.back().coordinates();
    lrange[0] = bodies.back().coordinates();

#pragma omp parallel
    {
      range_t trange;
      trange[1] = bodies.back().coordinates();
      trange[0] = bodies.back().coordinates();

#pragma omp parallel for
      for (size_t i = 0; i < bodies.size(); ++i) {
        for (size_t d = 0; d < gdimension; ++d) {
          if (bodies[i].coordinates()[d] + bodies[i].radius() > trange[1][d])
            trange[1][d] = bodies[i].coordinates()[d] + bodies[i].radius();
          if (bodies[i].coordinates()[d] - bodies[i].radius() < trange[0][d])
            trange[0][d] = bodies[i].coordinates()[d] - bodies[i].radius();
        } // for
      } // for
#pragma omp critical
      for (size_t d = 0; d < gdimension; ++d) {
        lrange[1][d] = std::max(lrange[1][d], trange[1][d]);
        lrange[0][d] = std::min(lrange[0][d], trange[0][d]);
      } // for 
    } // omp parallel 

    double max[gdimension];
    double min[gdimension];
    for (size_t i = 0; i < gdimension; ++i) {
      max[i] = lrange[1][i];
      min[i] = lrange[0][i];
    } // for

    MPI_Allreduce(MPI_IN_PLACE, max, dimension, MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, min, dimension, MPI_DOUBLE, MPI_MIN,
                  MPI_COMM_WORLD);

    for (size_t d = 0; d < gdimension; ++d) {
      range[0][d] = min[d];
      range[1][d] = max[d];
    } // for
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
      std::vector<body> &rbodies, const int64_t totalnbodies) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

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

    if (rank == 0) {
      master_recvcounts.resize(size);
    } // if

    // Echange the number of samples
    MPI_Gather(&nsample, 1, MPI_INT, &master_recvcounts[0], 1, MPI_INT, 0,
               MPI_COMM_WORLD);

    // Master
    // Sort the received keys and create the pivots
    if (rank == 0) {
      master_offsets.resize(size);
      master_nkeys = std::accumulate(master_recvcounts.begin(),
                                     master_recvcounts.end(), 0);
      if (totalnbodies < master_nkeys) {
        master_nkeys = totalnbodies;
      }
      // Number to receiv from each process
      for (int i = 0; i < size; ++i) {
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
    splitters.resize(size - 1 + 2);
    if (rank == 0) {
      std::sort(master_keys.begin(), master_keys.end(),
        [](auto &left, auto &right) {
          if (left.first < right.first) {
            return true;
          }
          if (left.first == right.first) {
            return left.second < right.second;
          }
          return false;
        }); // sort

      splitters[0].first = key_type::min();
      splitters[0].second = 0L;
      splitters[size].first = key_type::max();
      splitters[size].second = LONG_MAX;

      for (int i = 0; i < size - 1; ++i) {
        int64_t position = (master_nkeys / size) * (i + 1);
        splitters[i + 1] = master_keys[position];
        assert(splitters[i + 1].first > splitters[0].first &&
               splitters[i + 1].first < splitters[size].first);
      } // for
    } // if

    // Bradcast the splitters
    MPI_Bcast(&splitters[0],
              (size - 1 + 2) * sizeof(std::pair<key_type, int64_t>), MPI_BYTE,
              0, MPI_COMM_WORLD);
  }

private:
  // Key track of the splitter to know first and last key
  std::vector<std::pair<key_type, int64_t>> splitters_;

}; // class tree_colorer

#endif // _mpisph_tree_colorer_h_

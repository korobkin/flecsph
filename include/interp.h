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
 * @file interp.h
 * @author Oleg Korobkin
 * @date March 2022
 * @brief Various interpolation tools
 */

#pragma once

namespace interp {

typedef enum interp_grid_type_enum {
  uniform_grid,
  nonuniform_grid
} interp_grid_type;

/*
 * Unspecialized interpolating functor template class
 *
 * Template parameters:
 * - D: dimension
 * - G: grid kind, 0:uniform, 1:non-uniform
 *
 */
template<int D, interp_grid_type G> 
struct interpolating_function_u {};


/*
 * One-dimensional interpolator on a uniform grid
 *
 * Usage example:
 *   double ys[5] = {0.0, 0.1, 0.3, 0.6, 1.0};
 *   interp::interpolating_function_1d f{0, 1, ys, sizeof(ys)/sizeof(double)};
 *   std::cout << "f(0.215) = " << f(0.215) << std::endl; 
 *
 */
typedef struct interpolating_function_u<1, uniform_grid> interpolating_function_1d;
template<>
struct interpolating_function_u<1, uniform_grid> {

  // default constructor
  interpolating_function_u<1,uniform_grid>() {
    fs_ = nullptr;
    N_ = 0;
  }

  // constructor for a function with a uniform grid
  interpolating_function_u<1,uniform_grid>(const double x1, const double x2,
      const double * const fs, const int N) {
    fs_ = nullptr;
    set_data(x1, x2, fs, N);
  }

  void
  set_data(const double x1, const double x2,
      const double * const fs, const int N) {
    N_= N;
    if (fs_ != nullptr) free(fs_);
    fs_ = (double*) malloc(N*sizeof(double));
    x1_ = x1;
    x2_ = x2;
    dx_ = (x2 - x1)/(double)(N - 1);
    for (int i=0; i<N; i++) {
        fs_[i] = fs[i];
    }
  }

  ~interpolating_function_u<1,uniform_grid>() {
      if (fs_ != nullptr) free(fs_);
  }

  double
  operator() (const double x) {
    int i = (int) ((x - x1_)/dx_);
    if (i < 0) return fs_[0];
    if (i > N_ - 2) return fs_[N_-1];
    double xl = x1_ + i*dx_, yl = fs_[i];
    double xr = xl  + dx_,   yr = fs_[i+1];
    return yl + (x - xl)/dx_*(yr - yl);
  }

private:
    int N_;
    double x1_, x2_, dx_;
    double *fs_;
}; // interpolating_function<1, uniform_grid>



} // end namespace interp

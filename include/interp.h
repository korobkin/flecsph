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

/*
 * Interpolating functor class
 *
 * Template parameters:
 * - D: dimension
 * - G: grid kind, 0:uniform, 1:non-uniform
 *
 */
struct interpolating_function_1d {

  // default constructor
  interpolating_function_1d() {
    xs_ = nullptr;
    ys_ = nullptr;
    N_ = 0;
  }

  // constructor for a function with a uniform grid
  interpolating_function_1d(const double x1, const double x2,
      const double * const ys, const int N) {
    xs_ = nullptr;
    ys_ = nullptr;
    set_data(x1, x2, ys, N);
  }

  void
  set_data(const double x1, const double x2,
      const double * const ys, const int N) {
    N_= N;
    if (xs_ != nullptr) free(xs_);
    if (ys_ != nullptr) free(ys_);
    xs_ = (double*) malloc(N*sizeof(double));
    ys_ = (double*) malloc(N*sizeof(double));
    xs_[0]   = x1_ = x1;
    xs_[N-1] = x2_ = x2;
    dx_ = (x2 - x1)/(double)(N - 1);
    for (int i=0; i<N; i++) {
        xs_[i] = x1 + dx_*i;
        ys_[i] = ys[i];
    }
  }

  ~interpolating_function_1d() {
      if (xs_ != nullptr) free(xs_);
      if (ys_ != nullptr) free(ys_);
  }

  double
  operator() (const double x) {
    int i = (int) ((x - x1_)/dx_);
    if (i < 0) return ys_[0];
    if (i > N_ - 2) return ys_[N_-1];
    double xl = xs_[i],   yl = ys_[i];
    double xr = xs_[i+1], yr = ys_[i+1];
    return yl + (x - xl)/(xr - xl)*(yr - yl);
  }

private:
    int N_;
    double x1_, x2_, dx_;
    double *xs_, *ys_;
}; // interpolating_function_1d


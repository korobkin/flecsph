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
 * - T: class
 * - D: dimension
 * - O: interpolation order
 * - G: grid kind, 0:uniform, 1:non-uniform
 *
 */
template<typename T, int D, interp_grid_type G, int O>
struct interpolating_function_u {};


/*
 * One-dimensional interpolator on a uniform grid
 *
 * Usage example:
 *   double ys[5] = {0.0, 0.1, 0.3, 0.6, 1.0};
 *   interp::linear_interpolator_1d f{0, 1, ys, sizeof(ys)/sizeof(double)};
 *   std::cout << "f(0.215) = " << f(0.215) << std::endl;
 *
 */
typedef struct
interpolating_function_u<double, 1, uniform_grid, 1> linear_interpolator_1d;
typedef struct
interpolating_function_u<double, 1, uniform_grid, 3> cubic_interpolator_1d;

template<typename T, int O> struct
interpolating_function_u<T, 1, uniform_grid, O> {

  // default constructor
  interpolating_function_u<T, 1, uniform_grid, O>() {
    fs_ = nullptr;
    N_ = 0;
  }

  // constructor for a function with a uniform grid
  interpolating_function_u<T, 1, uniform_grid, O>(const double x1,
      const double x2, const T * const fs, const int N) {
    fs_ = nullptr;
    set_data(x1, x2, fs, N);
  }

  void
  set_data(const double x1, const double x2,
      const T * const fs, const int N) {
    N_= N;
    if (fs_ != nullptr) free(fs_);
    fs_ = (T*) malloc(N*sizeof(T));
    x1_ = x1;
    x2_ = x2;
    dx_ = (x2 - x1)/(double)(N - 1);
    for (int i=0; i<N; i++) {
        fs_[i] = fs[i];
    }
  }

  ~interpolating_function_u<T, 1, uniform_grid, O>() {
      if (fs_ != nullptr) free(fs_);
  }

  T operator() (const double x) {
    int i = (int) ((x - x1_)/dx_);
    T f{0};
    if (i < 0)
        f = fs_[0];
    else if (i > N_ - 2)
        f = fs_[N_-1];
    else {
        double xl = x1_ + i*dx_;
        double xr = xl  + dx_;
        T yl = fs_[i];
        T yr = fs_[i+1];
        if constexpr (O == 1) {
            assert(N_ > 1);
            f = yl + (x - xl)/dx_*(yr - yl);
        }

        if constexpr (O == 3) {
            assert(N_ > 3);
            if (i > 0) --i;
            if (i == N_- 3) --i;
            double xx0 = x - x1_ - i*dx_,
                   xx1 = xx0 - dx_,
                   xx2 = xx1 - dx_,
                   xx3 = xx2 - dx_;
            f = (fs_[i+3]*xx0 - fs_[i  ]*xx3)*xx1*xx2/3.
              - (fs_[i+2]*xx1 - fs_[i+1]*xx2)*xx0*xx3;
            f /= 2.*dx_*dx_*dx_;
        }
    }
    return f;
  }

private:
    int N_;
    double x1_, x2_, dx_;
    T *fs_;
}; // interpolating_function_u<T, 1, uniform_grid, O>

/**
* @brief   Finds an index i such that v[i]<= x < v[i+1]
*
* Returns: index i s.t. v[i] <= x < v[i+1]
*          if x < v[i], returns -1
*          if x >= v[N-1], returns N
* Note: vector v must be sorted, i.e. v[j]<= v[j+1]
*/
template<typename T> size_t
get_index(const T & x, const T * v, size_t v_size) {
  size_t i, i1, i2;
  i1 = 0;
  i2 = (v_size > 1) ? (v_size - 1) : 0;

  if (x*(1 + 1e-15) < v[0])
    i1 = -1;
  else if (x*(1 - 1e-15) > v[i2])
    i1 = i2;
  else {
    do{
      i = (i1 + i2)/2;
      T y = v[i];
      if (x < y)
        i2 = i;
      else
        i1 = i;
    } while(i2-i1>1);
  }
  return i1;
}

/**
* @brief   Finds an index i such that v[i]<= x < v[i+1]
*
* NOTE: Overload with the vector argument
*/
template<typename T> size_t
get_index(const T & x, const std::vector<T> & v) {
  return get_index(x, v.data(), v.size());
}


/*
 * One-dimensional interpolator on a nonuniform grid
 *
 * Usage example:
 *   double xs[] = {0.0, 0.5, 0.8, 0.9, 1.0};
 *   double ys[] = {0.0, 0.1, 0.3, 0.6, 1.0};
 *   interp::linear_interpolator_1d_nug f{xs, ys, sizeof(ys)/sizeof(double)};
 *   std::cout << "f(0.215) = " << f(0.215) << std::endl;
 *
 */
typedef struct
interpolating_function_u<double, 1, nonuniform_grid, 1> linear_interpolator_1d_nug;
typedef struct
interpolating_function_u<double, 1, nonuniform_grid, 3> cubic_interpolator_1d_nug;

template<typename T, int O>
struct interpolating_function_u<T, 1, nonuniform_grid, O> {

  // default constructor
  interpolating_function_u<T, 1, nonuniform_grid, O>() {
    xs_ = nullptr;
    fs_ = nullptr;
    N_ = 0;
  }

  interpolating_function_u<T, 1, nonuniform_grid, O>(const double * xs,
      const T * const fs, const int N) {
    xs_ = nullptr;
    fs_ = nullptr;
    set_data(xs, fs, N);
  }

  void
  set_data(const double * const xs, const T * const fs, const int N) {
    N_= N;
    if (xs_ != nullptr) free(xs_);
    if (fs_ != nullptr) free(fs_);
    xs_ = (double*) malloc(N*sizeof(double));
    fs_ = (T*) malloc(N*sizeof(T));
    for (int i=0; i<N; i++) {
        xs_[i] = xs[i];
        fs_[i] = fs[i];
    }
  }

  ~interpolating_function_u<T, 1, nonuniform_grid, O>() {
      if (xs_ != nullptr) free(xs_);
      if (fs_ != nullptr) free(fs_);
  }

  T operator() (const double x) {
    int i = get_index(x, xs_, N_);
    T f{0};
    if (i < 0)
        f = fs_[0];
    else if (i > N_ - 2)
        f = fs_[N_-1];
    else {
        if constexpr (O == 1) {
            assert(N_ > 1);
            int i1 = i, i2 = i + 1;
            double xx1 = x - xs_[i1],
                   xx2 = x - xs_[i2];
            f = (fs_[i1]*xx2 - fs_[i2]*xx1)/(xs_[i1] - xs_[i2]);
        }
        if constexpr (O == 3) {
            assert(N_ > 3);
            int i2 = i;
            if (i == 0)
                i2 = 1;
            if (i == N_ - 2)
                i2 = N_ - 3;

            int i1 = i2 - 1, i3 = i2 + 1, i4 = i2 + 2;

            double xx1 = x - xs_[i1],
                   xx2 = x - xs_[i2],
                   xx3 = x - xs_[i3],
                   xx4 = x - xs_[i4];

            f = (fs_[i1]/(xx2 - xx1)*xx4/(xx3 - xx1) - fs_[i4]/(xx2 - xx4)*xx1/(xx3 - xx4))
              * xx2*xx3/(xx4 - xx1)
              + (fs_[i2]/(xx1 - xx2)*xx3/(xx4 - xx2) - fs_[i3]/(xx1 - xx3)*xx2/(xx4 - xx3))
              * xx1*xx4/(xx3 - xx2);
        }
    }
    return f;
  }

private:
    int N_;
    double *xs_;
    T *fs_;
}; // interpolating_function_u<T, 1, nonuniform_grid, O> {


/*
 * Two-dimensional interpolator on a uniform grid
 *
 * Usage example:
 *   double fs[] = {0.0, 0.1, 0.3, // column-major order
 *                  0.6, 1.0, 2.0};
 *   interp::linear_interpolator_2d f{-1.0, 1.0, 0.1, 1.0, fs, 3, 2};
 *   std::cout << "f(-0.1, 0.215) = " << f(-0.1, 0.215) << std::endl;
 *
 */
typedef struct
interpolating_function_u<double, 2, uniform_grid, 1> linear_interpolator_2d;
typedef struct
interpolating_function_u<double, 2, uniform_grid, 3> cubic_interpolator_2d;

template<typename T, int O> struct
interpolating_function_u<T, 2, uniform_grid, O> {

  // default constructor
  interpolating_function_u<T, 2, uniform_grid, O>() {
    fs_ = nullptr;
    Nx_ = Ny_ = 0;
  }

  // constructor for a function with a uniform grid
  interpolating_function_u<T, 2, uniform_grid, O>(
      const double x1, const double x2,
      const double y1, const double y2,
      const T * const fs,
      const int Nx, const int Ny) {
    fs_ = nullptr;
    set_data(x1, x2, y1, y2, fs, Nx, Ny);
  }

  void
  set_data(
      const double x1, const double x2,
      const double y1, const double y2,
      const T * const fs,
      const int Nx, const int Ny) {
    Nx_= Nx; Ny_= Ny;
    x1_ = x1;  x2_ = x2;
    y1_ = y1;  y2_ = y2;
    if (fs_ != nullptr) free(fs_);
    fs_ = (T*) malloc(Nx*Ny*sizeof(T));
    dx_ = (x2 - x1)/(double)(Nx - 1);
    dy_ = (y2 - y1)/(double)(Ny - 1);
    for (int i=0; i<Nx*Ny; i++) {
        fs_[i] = fs[i];
    }
  }

  ~interpolating_function_u<T, 2, uniform_grid, O>() {
      if (fs_ != nullptr) free(fs_);
  }

  T operator() (const double x, const double y) {
    int i = (int) ((x - x1_)/dx_);
    i = (i == -1     && x*(1 + 1e-14) > x1_) ? 0  : i;
    i = (i == Nx_- 1 && x*(1 - 1e-14) < x2_) ? Nx_- 2 : i;

    int j = (int) ((y - y1_)/dy_);
    j = (j == -1     && y*(1 + 1e-14) > y1_) ? 0 : j;
    j = (j == Ny_- 1 && y*(1 - 1e-14) < y2_) ? Ny_- 2 : j;

    T f{0};
    if (0 <= i && i < Nx_- 1 && 0 <= j && j < Ny_- 1) {
        assert(Nx_ > 1 && Ny_ > 1);

        double h2 = (x - x1_)/dx_ - i, h1 = 1. - h2,
               g2 = (y - y1_)/dy_ - j, g1 = 1. - g2;

        int   i11 = i + j*Nx_, i12 = i11 + Nx_,
              i21 = i11 + 1,   i22 = i12 + 1;

        if constexpr (O == 1) {
            f = h1*(fs_[i11]*g1 + fs_[i12]*g2)
              + h2*(fs_[i21]*g1 + fs_[i22]*g2);
        }

        if constexpr (O == 3) {
            assert(Nx_ > 3 && Ny_ > 3);
            // TODO
        }
    }
    return f;
  }

private:
    int Nx_, Ny_;
    double x1_, x2_, dx_;
    double y1_, y2_, dy_;
    T *fs_;
}; // interpolating_function_u<T, 2, uniform_grid, O>


/*
 * Two-dimensional interpolator on a non-uniform grid
 *
 * Usage example:
 *   double xs[] = {-1.0, 0.5, 1.0}; // x-grid
 *   double ys[] = { 0.0, 0.1, 0.9, 1.0}; // y-grid
 *   double fs[] = {0.0, 0.1, 0.3, // column-major order
 *                  0.1,-1.0, 2.1,
 *                  0.3,-1.0, 2.0,
 *                  0.6, 0.3, 2.0};
 *   interp::linear_interpolator_2d_nug f{xs, ys, fs, 3, 4};
 *   std::cout << "f(-0.1, 0.215) = " << f(-0.1, 0.215) << std::endl;
 *
 */
typedef struct
interpolating_function_u<double, 2, nonuniform_grid, 1> linear_interpolator_2d_nug;
typedef struct
interpolating_function_u<double, 2, nonuniform_grid, 3> cubic_interpolator_2d_nug;

template<typename T, int O> struct
interpolating_function_u<T, 2, nonuniform_grid, O> {

  // default constructor
  interpolating_function_u<T, 2, nonuniform_grid, O>() {
    xs_ = nullptr;
    ys_ = nullptr;
    fs_ = nullptr;
    Nx_ = Ny_ = 0;
  }

  // constructor for a function with a uniform grid
  interpolating_function_u<T, 2, nonuniform_grid, O>(
      const double * const xs,
      const double * const ys,
      const T * const fs,
      const int Nx, const int Ny) {
    xs_ = nullptr;
    ys_ = nullptr;
    fs_ = nullptr;
    set_data(xs, ys, fs, Nx, Ny);
  }

  void
  set_data(
      const double * const xs,
      const double * const ys,
      const T * const fs,
      const int Nx, const int Ny) {
    Nx_= Nx; Ny_= Ny;
    if (xs_ != nullptr) free(xs_);
    if (ys_ != nullptr) free(ys_);
    if (fs_ != nullptr) free(fs_);
    xs_ = (double*) malloc(Nx * sizeof(double));
    ys_ = (double*) malloc(Ny * sizeof(double));
    fs_ = (T*) malloc(Nx * Ny * sizeof(T));
    memcpy(xs_, xs, Nx * sizeof(double));
    memcpy(ys_, ys, Ny * sizeof(double));
    memcpy(fs_, fs, Nx * Ny * sizeof(T));
  }

  ~interpolating_function_u<T, 2, nonuniform_grid, O>() {
      if (xs_ != nullptr) free(xs_);
      if (ys_ != nullptr) free(ys_);
      if (fs_ != nullptr) free(fs_);
  }

  T operator() (const double x, const double y) {
    int i = get_index(x, xs_, Nx_);
    int j = get_index(y, ys_, Ny_);

    T f{0};
    if (0 <= i && i < Nx_- 1 && 0 <= j && j < Ny_- 1) {

        double h2 = (x - xs_[i])/(xs_[i+1] - xs_[i]),
               h1 = 1. - h2,
               g2 = (y - ys_[j])/(ys_[j+1] - ys_[j]),
               g1 = 1. - g2;

        int   i11 = i + j*Nx_, i12 = i11 + Nx_,
              i21 = i11 + 1,   i22 = i12 + 1;

        if constexpr (O == 1) {
            assert(Nx_ > 1 && Ny_ > 1);
            f = h1*(fs_[i11]*g1 + fs_[i12]*g2)
              + h2*(fs_[i21]*g1 + fs_[i22]*g2);
        }

        if constexpr (O == 3) {
            assert(Nx_ > 3 && Ny_ > 3);
            // TODO
        }
    }
    return f;
  }

private:
    int Nx_, Ny_;
    double *xs_, *ys_;
    T *fs_;
}; // interpolating_function_u<T, 2, nonuniform_grid, O>


/*
 * Three-dimensional interpolator on a uniform grid
 *
 * Usage example:
 * double fs[] = {0.1, -0.1, -0.3,
 *                0.6,  1.0,  2.4,
 *                0.3,  1.3,  2.1,
 *
 *               -0.4, -0.8,  1.3,
 *                0.6,  0.2,  0.9,
 *                0.3,  0.3,  0.3};
 * interp::linear_interpolator_3d f{-1., 1., 0., 1., 0., 2.5,  fs, 3, 3, 2};
 * 
 */
typedef struct
interpolating_function_u<double, 3, uniform_grid, 1> linear_interpolator_3d;
typedef struct
interpolating_function_u<double, 3, uniform_grid, 3> cubic_interpolator_3d;

template<typename T, int O> struct
interpolating_function_u<T, 3, uniform_grid, O> {

  // default constructor
  interpolating_function_u<T, 3, uniform_grid, O>() {
    fs_ = nullptr;
    Nx_ = Ny_ = Nz_ = 0;
  }

  // constructor for a function with a uniform grid
  interpolating_function_u<T, 3, uniform_grid, O>(
      const double x1, const double x2,
      const double y1, const double y2,
      const double z1, const double z2,
      const T * const fs,
      const int Nx, const int Ny, const int Nz) {
    fs_ = nullptr;
    set_data(x1, x2, y1, y2, z1, z2, fs, Nx, Ny, Nz);
  }

  void
  set_data(
      const double x1, const double x2,
      const double y1, const double y2,
      const double z1, const double z2,
      const T * const fs,
      const int Nx, const int Ny, const int Nz) {
    Nx_= Nx; Ny_= Ny; Nz_ = Nz;
    x1_ = x1;  x2_ = x2;
    y1_ = y1;  y2_ = y2;
    z1_ = z1;  z2_ = z2;
    if (fs_ != nullptr) free(fs_);
    fs_ = (T*) malloc(Nx*Ny*Nz*sizeof(T));
    dx_ = (x2 - x1)/(double)(Nx - 1);
    dy_ = (y2 - y1)/(double)(Ny - 1);
    dz_ = (z2 - z1)/(double)(Nz - 1);
    for (int i=0; i<Nx*Ny*Nz; i++) {
        fs_[i] = fs[i];
    }
  }

  ~interpolating_function_u<T, 3, uniform_grid, O>() {
      if (fs_ != nullptr) free(fs_);
  }

  T operator() (const double x, const double y, const double z) {
    int i = (int) ((x - x1_)/dx_);
    i = (i == -1     && x*(1 + 1e-14) > x1_) ? 0  : i;
    i = (i == Nx_- 1 && x*(1 - 1e-14) < x2_) ? Nx_- 2 : i;

    int j = (int) ((y - y1_)/dy_);
    j = (j == -1     && y*(1 + 1e-14) > y1_) ? 0 : j;
    j = (j == Ny_- 1 && y*(1 - 1e-14) < y2_) ? Ny_- 2 : j;

    int k = (int) ((z - z1_)/dz_);
    k = (k == -1     && z*(1 + 1e-14) > z1_) ? 0 : k;
    k = (k == Nz_- 1 && z*(1 - 1e-14) < z2_) ? Nz_- 2 : k;

    T f{0};
    if (0 <= i && i < Nx_- 1 && 0 <= j && j < Ny_- 1 && 0 <= k && k < Nz_-1) {
        assert(Nx_ > 1 && Ny_ > 1 && Nz_ > 1);

        double h2 = (x - x1_)/dx_ - i, h1 = 1. - h2,
               g2 = (y - y1_)/dy_ - j, g1 = 1. - g2,
	       m2 = (z - z1_)/dz_ - k, m1 = 1. - m2;

	int   i111 = i + j*Ny_ + k*Ny_*Nz_,
	      i121 = i111 + Nx_,
	      i112 = i111 + Nx_*Ny_,
	      i122 = i121 + Nx_*Ny_,
	      i211 = i111 + 1,
	      i221 = i121 + 1,
	      i212 = i112 + 1,
	      i222 = i122 + 1;

	if constexpr (O == 1) {
	    f = m1*((fs_[i111]*h1 + fs_[i211]*h2)*g1 
	          + (fs_[i121]*h1 + fs_[i221]*h2)*g2) 
	      + m2*((fs_[i112]*h1 + fs_[i212]*h2)*g1 
	          + (fs_[i122]*h1 + fs_[i222]*h2)*g2);
	}	

        if constexpr (O == 3) {
            assert(Nx_ > 3 && Ny_ > 3);
            // TODO
        }
    }
    return f;
  }

private:
    int Nx_, Ny_, Nz_;
    double x1_, x2_, dx_;
    double y1_, y2_, dy_;
    double z1_, z2_, dz_;
    T *fs_;
}; // interpolating_function_u<T, 3, uniform_grid, O>


/*
 * Three-dimensional interpolator on a non-uniform grid
 *
 * Usage example:
 *   double xs[] = {-1.0, 0.5, 1.0}; // x-grid
 *   double ys[] = { 0.0, 0.1, 0.9, 1.0}; // y-grid
 *   double zs[] = { 0.0, 1.0}; // z-grid
 *   double fs[] = {0.0, 0.1, 0.3, // column-major order
 *                  0.1,-1.0, 2.1,
 *                  0.3,-1.0, 2.0,
 *                  0.6, 0.3, 2.0};
 *   interp::linear_interpolator_3d_nug f{xs, ys, zs, fs, 3, 4, 2};
 *   std::cout << "f(-0.1, 0.215) = " << f(-0.1, 0.215) << std::endl;
 *
 */
typedef struct
interpolating_function_u<double, 3, nonuniform_grid, 1> linear_interpolator_3d_nug;
typedef struct
interpolating_function_u<double, 3, nonuniform_grid, 3> cubic_interpolator_3d_nug;

template<typename T, int O> struct
interpolating_function_u<T, 3, nonuniform_grid, O> {

  // default constructor
  interpolating_function_u<T, 3, nonuniform_grid, O>() {
    xs_ = nullptr;
    ys_ = nullptr;
    zs_ = nullptr;
    fs_ = nullptr;
    Nx_ = Ny_ = Nz_ = 0;
  }

  // constructor for a function with a uniform grid
  interpolating_function_u<T, 3, nonuniform_grid, O>(
      const double * const xs,
      const double * const ys,
      const double * const zs,
      const T * const fs,
      const int Nx, const int Ny, const int Nz) {
    xs_ = nullptr;
    ys_ = nullptr;
    zs_ = nullptr;
    fs_ = nullptr;
    set_data(xs, ys, zs, fs, Nx, Ny, Nz);
  }

  void
  set_data(
      const double * const xs,
      const double * const ys,
      const double * const zs,
      const T * const fs,
      const int Nx, const int Ny, const int Nz) {
    Nx_= Nx; Ny_= Ny; Nz_ = Nz;
    if (xs_ != nullptr) free(xs_);
    if (ys_ != nullptr) free(ys_);
    if (zs_ != nullptr) free(zs_);
    if (fs_ != nullptr) free(fs_);
    xs_ = (double*) malloc(Nx * sizeof(double));
    ys_ = (double*) malloc(Ny * sizeof(double));
    zs_ = (double*) malloc(Nz * sizeof(double));
    fs_ = (T*) malloc(Nx * Ny * Nz * sizeof(T));
    memcpy(xs_, xs, Nx * sizeof(double));
    memcpy(ys_, ys, Ny * sizeof(double));
    memcpy(zs_, zs, Nz * sizeof(double));
    memcpy(fs_, fs, Nx * Ny * sizeof(T));
  }

  ~interpolating_function_u<T, 3, nonuniform_grid, O>() {
      if (xs_ != nullptr) free(xs_);
      if (ys_ != nullptr) free(ys_);
      if (zs_ != nullptr) free(zs_);
      if (fs_ != nullptr) free(fs_);
  }

  T operator() (const double x, const double y, const double z) {
    int i = get_index(x, xs_, Nx_);
    int j = get_index(y, ys_, Ny_);
    int k = get_index(z, zs_, Nz_);
   
    //printf("%d %d %d\n", i, j, k);
    T f{0};
    if (0 <= i && i < Nx_- 1 && 0 <= j && j < Ny_- 1 && 0 <= k && k < Nz_-1) {
        assert(Nx_ > 1 && Ny_ > 1 && Nz_ > 1);

        double h2 = (x - xs_[i])/(xs_[i+1] - xs_[i]), h1 = 1. - h2,
               g2 = (y - ys_[j])/(ys_[j+1] - ys_[j]), g1 = 1. - g2,
	       m2 = (z - zs_[k])/(zs_[k+1] - zs_[k]), m1 = 1. - m2;

	int   i111 = i + j*Ny_ + k*Ny_*Nz_,
	      i121 = i111 + Nx_,
	      i112 = i111 + Nx_*Ny_,
	      i122 = i121 + Nx_*Ny_,
	      i211 = i111 + 1,
	      i221 = i121 + 1,
	      i212 = i112 + 1,
	      i222 = i122 + 1;

	if constexpr (O == 1) {
	    f = m1*((fs_[i111]*h1 + fs_[i211]*h2)*g1 
	          + (fs_[i121]*h1 + fs_[i221]*h2)*g2) 
	      + m2*((fs_[i112]*h1 + fs_[i212]*h2)*g1 
	          + (fs_[i122]*h1 + fs_[i222]*h2)*g2);
	}	

        if constexpr (O == 3) {
            assert(Nx_ > 3 && Ny_ > 3);
            // TODO
        }
    }
    return f;
  }

private:
    int Nx_, Ny_, Nz_;
    double *xs_, *ys_, *zs_;
    T *fs_;
}; // interpolating_function_u<T, 3, nonuniform_grid, O>
} // end namespace interp

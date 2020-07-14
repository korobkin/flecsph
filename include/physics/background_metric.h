/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2020 Triad National Security, LLC
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
 * @file background_metric.h
 * @brief Define several different background metric
 *        to compute geneneral relativistic accleration.
 *        We follow (-1,1,1,1) signature
 */

#pragma once

#include "params.h"
#include "tensor.h"

//NOTE : Define this as namespace??

double gc = param::gravitational_constant;

using sym_tensor_rank2 = flecsi::tensor_u<double, symmetry_type::symmetric, 4, 4>;
using sym_tensor_rank3 = flecsi::tensor_u<double, symmetry_type::symmetric, 4, 4, 4>;

// Flat Minkowski metric in Carteisan coordinate
// Rank = 2, Dim = 4 -> 10 independent quantities
sym_tensor_rank2 gMinkowski{0};

gMinkowski(0,0) = -1.0; //tt
gMinkowski(0,0) =  0.0; //tx
gMinkowski(0,2) =  0.0; //ty
gMinkowski(0,3) =  0.0; //tz
gMinkowski(1,1) =  1.0; //xx
gMinkowski(1,2) =  0.0; //xy
gMinkowski(1,3) =  0.0; //xz
gMinkowski(2,2) =  1.0; //yy
gMinkowski(2,3) =  0.0; //yz
gMinkowski(3,3) =  1.0; //zz

// First derivative of metric
// It is used to compute acceleration equation
// Rank = 3, Dim = 4 -> 20 independent quantities
sym_tensor_rank3 d_gMinkowski{0};

d_gMinkowski(0,0,0) = 0.0; //ttt
d_gMinkowski(0,0,1) = 0.0; //ttx
d_gMinkowski(0,0,2) = 0.0; //tty
d_gMinkowski(0,0,3) = 0.0; //ttz
d_gMinkowski(0,1,1) = 0.0; //txx
d_gMinkowski(0,1,2) = 0.0; //txy
d_gMinkowski(0,1,3) = 0.0; //txz
d_gMinkowski(0,2,2) = 0.0; //tyy
d_gMinkowski(0,2,3) = 0.0; //tyz
d_gMinkowski(0,3,3) = 0.0; //tzz
d_gMinkowski(1,1,1) = 0.0; //xxx
d_gMinkowski(1,1,2) = 0.0; //xxy
d_gMinkowski(1,1,3) = 0.0; //xxz
d_gMinkowski(1,2,2) = 0.0; //xyy
d_gMinkowski(1,2,3) = 0.0; //xyz
d_gMinkowski(1,3,3) = 0.0; //xzz
d_gMinkowski(2,2,2) = 0.0; //yyy
d_gMinkowski(2,2,3) = 0.0; //yyz
d_gMinkowski(2,3,3) = 0.0; //yzz
d_gMinkowski(3,3,3) = 0.0; //zzz

// Static spherically symmetric metric (i.e. Schwarzschild) in Kerr-Schild coordinate
sym_tensor_rank2 gSchwarz{0};

// Define coordinates and physical quantities
//NOTE : How we understand this quantities? 
//       This shouldn't be realted with particles' evolution
//TODO : Change it to relevant form. Save it as now to get clear view

const double M_back = 1.0;
double coords[4] = {0}; //General spacetime coordiantes
double x = coords[0], y = coords[1], z = coords[2]; // short hand notation for spatial coordinates
double r = std::sqrt(x*x + y*y + z*z);
//It is good to define lapse and shift to simplify expression

double alpha = std::sqrt(r/(r+2*M_back));

double beta_u[3], beta_d[3], beta_sum;

for(int i = 1; i < 4; ++i) {
  beta_u[i] = 2*M_back/r*coords[i]/(r+2*M_back);
  beta_d[i] = 2*M_back*coords[i]/(r*r);
  beta_sum += beta_u[i]*beta_d[i];
}

gSchwarz(0,0) = -alpha*alpha + beta_sum; //tt

//tx, ty, tz components
for(int i = 1; i < 4; ++i) {
  gSchwarz(0,i) = beta_u[i];
}

// ij (spatial) components
for(int i = 1; i < 4; ++i) {
  for(int j = i; j < 4; ++j) {
    gSchwarz(i,j) = (i==j) + 2*M_back*coords[i]*coords[j]/(r*r*r);
  }
}


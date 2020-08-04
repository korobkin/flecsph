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
using gen_tensor_rank2 = flecsi::tensor_u<double, symmetry_type::generic, 4, 4>;
using gen_tensor_rank3 = flecsi::tensor_u<double, symmetry_type::generic, 4, 4, 4>;

//HL : I am not using loop here to see indices. Will simplify later

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

constexpr double C_LIGHT_CGS = 2.99792458e10; // Speed of light in CGS
const double M_back = 1.0;
const double r_sch = 2.*M_back*gc/(C_LIGHT_CGS*C_LIGHT_CGS); // Schwarzschild radius
double coords[4] = {0}; //General spacetime coordiantes
double x = coords[0], y = coords[1], z = coords[2]; // short hand notation for spatial coordinates
double r = std::sqrt(x*x + y*y + z*z);
double r2 = r*r;
double r3 = r*r2;
double r4 = r2*r2;
double r5 = r2*r3;
//It is good to define lapse and shift to simplify expression

double alpha = std::sqrt(r/(r+r_sch));

double beta_u[4], beta_d[4], beta_sum;

beta_u[0] = 0.0, beta_d[0] = 0.0;

for(int i = 1; i < 4; ++i) {
  beta_u[i] = r_sch/r*coords[i]/(r+2*M_back);
  beta_d[i] = r_sch*coords[i]/(r2);
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
    gSchwarz(i,j) = (i==j) + r_sch*coords[i]*coords[j]/(r3);
  }
}

//First derivative of metric. Now this is not fully symmetric
gen_tensor_rank3 d_gSchwarz{0};

//tab components. All zeros
for(int i = 0; i < 4; ++i) {
  for(int j = 0; j < 4; ++j) {
    d_gSchwarz(0,i,j) = 0.0;
  }
}

//itj (or ijt) components 
for(int i = 1; i < 4; ++i) {
  for(int j = 1; j < 4; ++j) {
    d_gSchwarz(i,0,j) = r_sch*(i==j)/(r*(r+2*M_back));
    d_gSchwarz(i,j,0) = d_gSchwarz(i,0,j);
  }
}

//ijk (all spatial) components
for(int i = 1; i < 4; ++i) {
  for(int j = 1; j < 4; ++j) {
    for(int k = 1; k < 4; ++k) {
      d_gSchwarz(i,j,k) = r_sch/(r3)*((i==k) + (j==k)) 
                           + 3*r_sch*coords[i]*coords[j]*coords[k]/r5;
    }
  }
}

// Axisymmetic metric (i.e. Kerr) in Kerr-Schild Cartesian coordinate

// NOTE : I keep both Schwarzschild and Kerr for now for sanity check. 
//        Once everything looks fine, I will remove Schwarzschild since
//        a->0 (or J->0) in Kerr will return Schwarzschild

sym_tensor_rank2 gKerr{0};

//Define dimensionaless spin
//TODO : make it as parameter
const double J_ang = 0.1; //Angular momentum
const double a_ang = J_ang/(M_back*C_LIGHT_CGS); // Spin parameter

//Define scalar quantities 

double f_scalar = 2.*gc*M_back*r3/(r4+a_ang*a_ang*z*z);

//Define k 4-vector in covariant form
double k_vec[4];
k_vec[0] = 1.0;
k_vec[1] = (r*x+a_ang*y)/(r2+a_ang*a_ang); 
k_vec[2] = (r*y-a_ang*x)/(r2+a_ang*a_ang); 
k_vec[3] = z/r; 

for(int i = 0; i < 4; ++i){
  for(int j = 0; j < 4; ++j){
    gKerr(i,j) = gMinkowski(i,j) + f_scalar*k_vec[i]*k_vec[j];
  }
}

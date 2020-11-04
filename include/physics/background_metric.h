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

//#include "params.h"
//#include "tensor.h"

namespace background_metric{
#include "tensor.h"
#include "params.h"
#include "phys_consts.h"
#include "density_profiles.h"



void reset_metric(sym_tensor_rank2_spacetime & gm,
                  sym_tensor_rank2_spacetime & inv_gm,
                  sym_tensor_rank2_spacetime (&d_gm)[4]){
  //log_one(info) << "reset the metric" << std::endl;
  for(int mu = 0; mu < 4; ++mu) {
    for(int nu = mu; nu < 4; ++nu) {
      gm(mu,nu) = 0;
      inv_gm(mu,nu) = 0;
      for(int l = 0; l < 4; ++l) {
        (d_gm[l])(mu,nu) = 0;
      }
    }
  }
}

void set_Minkowski_metric(const point_t & pos,
                          sym_tensor_rank2_spacetime & gMin,
                          sym_tensor_rank2_spacetime & inv_gMin,
                          sym_tensor_rank2_spacetime (&d_gMin)[4]){
  reset_metric(gMin, inv_gMin, d_gMin);
  gMin(0,0) += -1;
  inv_gMin(0,0) += -1;
  for(int i = 1; i < 4; ++i) {
    gMin(i,i) = 1;
    inv_gMin(i,i) = 1;
  }
  //d_gMin[1](1,1) = 1.0;
  /*for(int i = 1; i < 4; i++){
    for(int j = i; j < 4; j++){
      for(int k = 1; k < 4; ++k) {
        // TODO very weird bug: d_gm cannot be 0 or the driver never advance
        (d_gMin[k])(i,j) = 0.0e-90;
      }
    }
  }*/
#if 0
  log_one(info)<<"The metric: "<<std::endl;
  for(int mu = 0; mu < 4; mu ++) {
    for(int nu = 0; nu < 4; nu ++){
      log_one(info)<< "g"<<mu<<nu<<" " << gMin(mu,nu) << std::endl;
    }
  }
  for(int mu = 0; mu < 4; mu ++) {
    for(int nu = 0; nu < 4; nu ++){
      log_one(info)<< "g_inv"<<mu<<nu<<" " << inv_gMin(mu,nu) << std::endl;
    }
  }
#endif 

}  

void set_TOV_metric(const point_t & pos,
                    sym_tensor_rank2_spacetime & gTOV,
                    sym_tensor_rank2_spacetime & inv_gTOV,
                    sym_tensor_rank2_spacetime (&d_gTOV)[4]){
  double alpha2, dalpha2dr, beta2, dbeta2dr;
  double x = pos[0], y = pos[1], z = pos[2]; // short hand notation for spatial coordinates
  double R_star = param::sphere_radius / L_GEOM_TO_CGS;
  double rmin = R_star * 1.0e-15;
  double r = std::sqrt(x*x + y*y + z*z) + rmin; // Radial distance that we will use
  double xi[4] = {0, x/r, y/r, z/r}; // dimensionless coordinates
  double r2 = r*r;
  reset_metric(gTOV, inv_gTOV, d_gTOV);
  
  // In the interior of the star: TOV metric
  if(r < R_star){
    alpha2 = density_profiles::spherical_alpha2(r / R_star);
    dalpha2dr = density_profiles::spherical_dalpha2_dr(r / R_star);
    beta2 = density_profiles::spherical_beta2(r / R_star);
    dbeta2dr = density_profiles::spherical_dbeta2_dr(r / R_star);
    //log_one(info) << "alpha: "<<alpha2 <<"beta:" << beta2<< std::endl;
  } else {
    double r_sch = 2*param::sphere_mass/M_GEOM_TO_CGS;
    alpha2 = (1 - r_sch / r);
    dalpha2dr = r_sch/r2;
    beta2 = r/ (r - r_sch);
    dbeta2dr = -r_sch / (r - r_sch);
    //log_one(info) << "alpha: "<<alpha2 <<"beta:" << beta2<< std::endl;
  }
  gTOV(0,0) = -alpha2;
  inv_gTOV(0,0) = -1.0/alpha2;
  for(int i = 1; i< 4; i++) {
    (d_gTOV[i])(0,0) = -dalpha2dr * xi[i];
  }
  for(int i = 1; i < 4; i++){
    for(int j = 1; j < 4; j++){
      gTOV(i,j)       = (int)(i==j) + xi[i] * xi[j] * (beta2-1);
      inv_gTOV(i,j)   = (int)(i==j) + xi[i] * xi[j] * (1.0/beta2-1);
      for(int k = 1; k < 4; ++k) {
        (d_gTOV[k])(i,j) = ((int)(i==k)*xi[j] + (int)(j==k)*xi[i] - 2*xi[i]*xi[j]*xi[k])
                           *(beta2 - 1)/r
                           + xi[i]*xi[j]*(dbeta2dr);
      }
    }
  }
} //set_TOV_metric

} //namespace background_metric

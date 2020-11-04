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
 *        We follow (-1,1,1,1) signature.
 *        All quantities are expressed in geometrized unit
 */

#pragma once

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

//TODO : Check BH backgroud metics into Sparta to see Ham and mom constraint values

// Static spherically symmetric metric (i.e. Schwarzschild) in Kerr-Schild coordinate
void set_Schwarzschild_metric_KS(const point_t &pos,
                                 sym_tensor_rank2_spacetime & gSchwarzKS,
                                 sym_tensor_rank2_spacetime (&d_gSchwarzKS)[4]){

//Some precomputation
constexpr double M_BH = 1.0; // Background BH metric mass
constexpr double r_sch = 2.*M_BH; // Schwarzschild radius in geometrized unit
double x = pos[0], y = pos[1], z = pos[2]; // short hand notation for spatial coordinates
double coords[4] = {0.0, x, y, z}; //General spacetime coordiantes
double r_floor = 1e-6; // Small floor value to avoid radial distance r goes to zero
double r = std::sqrt(x*x + y*y + z*z) + r_floor; // Radial distance that we will use TODO : Maybe not a great idea...
double r2 = r*r;
double r3 = r*r2;
double r4 = r2*r2;
double r5 = r2*r3;

//It is good to define lapse and shift to simplify expression
double alpha = std::sqrt(r/(r+r_sch));
double beta_u[4], beta_d[4], beta_sum;

beta_u[0] = 0.0, beta_d[0] = 0.0;

for(int i = 1; i < 4; ++i) {
  beta_u[i] = r_sch/r*coords[i]/(r+r_sch);
  beta_d[i] = r_sch*coords[i]/(r2);
  beta_sum += beta_u[i]*beta_d[i];
}

gSchwarzKS(0,0) = -alpha*alpha + beta_sum; //tt

//tx, ty, tz components
for(int i = 1; i < 4; ++i) {
  gSchwarzKS(0,i) = beta_d[i];
}

// ij (spatial) components
for(int i = 1; i < 4; ++i) {
  for(int j = 1; j < 4; ++j) {
    gSchwarzKS(i,j) = (i==j) + r_sch*coords[i]*coords[j]/(r3);
  }
}

//First derivative of metric. 

// tab components. All zeros
for(int i = 0; i < 4; ++i) {
  for(int j = 0; j < 4; ++j) {
    d_gSchwarzKS[0](i,j) = 0.0;
  }
}

// itt components
for(int i = 1; i < 4; ++i){
  d_gSchwarzKS[i](0,0) = coords[i]/((r+r_sch)*(r+r_sch)) - coords[i]/(r*(r+r_sch))
                       - r_sch*r_sch*coords[i]*(2.0*r+r_sch)/(r3*(r+r_sch)*(r+r_sch));
}

// itj (or ijt) components
for(int i = 1; i < 4; ++i) {
  for(int j = 1; j < 4; ++j) {
    d_gSchwarzKS[i](0,j) = r_sch*(i==j)/(r*(r+r_sch));
  }
}

// ijk (all spatial) components
for(int i = 1; i < 4; ++i) {
  for(int j = 1; j < 4; ++j) {
    for(int k = 1; k < 4; ++k) {
      d_gSchwarzKS[i](j,k) = r_sch/(r3)*((i==k) + (j==k)) 
                           + 3*r_sch*coords[1]*coords[j]*coords[k]/r5;
    }
  }
}

} // Schwarzschild in KS

// Schwarzschild in isotropic Cartesian coordinates
void set_Schwarzschild_metric_BL(const point_t &pos,
                                 sym_tensor_rank2_spacetime & gSchwarzIso,
                                 sym_tensor_rank2_spacetime (&d_gSchwarzIso)[4]){

constexpr double M_BH = 1.0; // Background BH metric mass
constexpr double r_sch = 2.*M_BH; // Schwarzschild radius in geometrized unit
double x = pos[0], y = pos[1], z = pos[2]; // short hand notation for spatial coordinates
double coords[4] = {0.0, x, y, z}; //General spacetime coordiantes
double r_floor = 1e-6; // Small floor value to avoid radial distance r goes to zero
double r = std::sqrt(x*x + y*y + z*z) + r_floor; // Radial distance that we will use TODO : Maybe not a great idea...
double r2 = r*r;
double r3 = r2*r;

// Prefactors
double f1 = 1.0 - r_sch/(4.0*r);
double f2 = 1.0 + r_sch/(4.0*r);

double f1sq = f1*f1, f2sq = f2*f2;
double f2qu = f2sq*f2sq;

// tt component
gSchwarzIso(0,0) = -f1sq/f2sq;

// ii (spatial diagonal) components
for(int i = 1; i < 4; ++i) {
  gSchwarzIso(i,i) = f2qu;
}

// Off diagonal components
for(int i = 0; i < 4; ++i){
  for(int j = i+1; j < 4; ++j){
    gSchwarzIso(i,j) = 0.0;
  }
}

// First derivative of metric

// tab component 
for(int i = 0; i < 4; ++i) {
  for(int j = 0; j < 4; ++j) {
    d_gSchwarzIso[0](i,j) = 0.0;
  }
}

// itt components
for(int i = 1; i < 4; ++i){
  d_gSchwarzIso[i](0,0) = -r_sch*coords[i]*f1sq/(2.0*r3*f2sq*f2) - r_sch*coords[i]*f1/(2.0*r3*f2sq);
}

// itj (or ijt) components
for(int i = 1; i < 4; ++i) {
  for(int j = 1; j < 4; ++j) {
    d_gSchwarzIso[i](0,j) = 0.0;
  }
}

// ijj components
for(int i = 1; i < 4; ++i) {
  for(int j = 1; j < 4; ++j) {
    d_gSchwarzIso[i](j,j) = - r_sch*coords[j]/r3*f2sq*f2 ;
  }
}

// i(spatial off-diagonal) components
for(int i = 1; i < 4; ++i){
  for(int j = 1; j < 4; ++j){
    for(int k = i+1; k < 4; ++k){
      d_gSchwarzIso[i](j,k) = 0.0;
    }
  }
}

} // Schwarzschild in Iso

// Axisymmetic metric (i.e. Kerr) in Kerr-Schild Cartesian coordinate
// NOTE : I keep both Schwarzschild and Kerr for now for sanity check. 
//        Once everything looks fine, I will remove Schwarzschild since
//        a->0 (or J->0) in Kerr will return Schwarzschild
void set_Kerr_metric_KS(const point_t &pos,
                        sym_tensor_rank2 & gKerrKS,
                        sym_tensor_rank2 (&d_gKerrKS)[4]){

constexpr double M_BH = 1.0; // Background BH metric mass
constexpr double r_sch = 2.*M_BH; // Schwarzschild radius in geometrized unit
double x = pos[0], y = pos[1], z = pos[2]; // short hand notation for spatial coordinates
double coords[4] = {0.0, x, y, z}; //General spacetime coordiantes
double r_floor = 1e-6; // Small floor value to avoid radial distance r goes to zero
double r = std::sqrt(x*x + y*y + z*z) + r_floor; // Radial distance that we will use TODO : Maybe not a great idea...
double r2 = r*r;
double r3 = r2*r;
double r4 = r2*r2;
double r5 = r2*r3;
// Define dimensionaless spin
// TODO : make it as parameter
const double J_ang = 0.1; //Angular momentum
const double a_ang = J_ang/M_BH; // Spin parameter in geometrical unit

// Some short hand notation
const double a2 = a_ang*a_ang;
double x2 = x*x, y2 = y*y, z2 = z*z;

// Define scalar quantities 
double f_scalar = 2.*GNEWT*r3/(r4+a2*z2); // in geometrical unit

// Define k 4-vector in covariant form
double k_vec[4];
k_vec[0] = 1.0;
k_vec[1] = (r*x+a_ang*y)/(r2+a2); 
k_vec[2] = (r*y-a_ang*x)/(r2+a2); 
k_vec[3] = z/r; 

// tt component
gKerrKS(0,0) = -1.0 + f_scalar;
// ij component
for(int i = 1; i < 4; ++i){
  for(int j = 1; j < 4; ++j){
    gKerrKS(i,j) = (i==j) + f_scalar*k_vec[i]*k_vec[j];
  }
}

// First derivative of metric

// Define derivative qunatities
// Derivatives of scalar
double d_f[4];
d_f[0] = 0.0;
d_f[1] = 6.0*M_BH*x*r/(a2*z2+r4)-8.0*M_BH*x*r5/((a2*z2+r4)*(a2*z2+r4));
d_f[2] = 6.0*M_BH*y*r/(a2*z2+r4)-8.0*M_BH*y*r5/((a2*z2+r4)*(a2*z2+r4));
d_f[3] = 6.0*M_BH*y*r/(a2*z2+r4)-2.0*M_BH*r3*(2.0*a2*z+4.0*z*r2)/((a2*z2+r4)*(a2*z2+r4));
// Derivatives of k-vector
double d_k[4][4];
// ta or at components are zero
for(int i = 1; i < 4; ++i){
   d_k[0][0] = 0.0;
   d_k[0][i] = 0.0;
   d_k[i][0] = 0.0;

}
// ij (spatial) components
d_k[1][1] = (x2/r+r)/(a2+r2)-2*x*(a_ang*y+x*r)/((a2+r2)*(a2+r2));
d_k[1][2] = (x*y/r+a_ang)/(a2+r2)-2*y*(a_ang*y+x*r)/((a2+r2)*(a2+r2));
d_k[1][3] = x*z/(r*(a2+r2))-2*z*(a_ang*y+x*r)/((a2+r2)*(a2+r2));
d_k[2][1] = (x*y/r-a_ang)/(a2+r2)-2*x*(y*r-a_ang*x)/((a2+r2)*(a2+r2));
d_k[2][2] = (y2/r+r)/(a2+r2)-2*y*(y*r-a_ang*x)/((a2+r2)*(a2+r2));
d_k[2][3] = y*z/(r*(a2+r2))-2*z*(y*r-a_ang*x)/((a2+r2)*(a2+r2));
d_k[3][1] = -x*z/r3;
d_k[3][2] = -y*z/r3;
d_k[3][3] = -z2/r3 + 1/r;

// Define first derivative of Kerr metric
for(int i = 0; i < 4; ++i) {
  for(int j = 0; j < 4; ++j) {
    for(int k = 0; k < 4; ++k) {
      d_gKerrKS[i](j,k) = d_f[i]*k_vec[j]*k_vec[k] 
                      + f_scalar*(d_k[k][i]*k_vec[j] + k_vec[i]*d_k[k][j]);
    }
  }
}

} // Kerr in KS


// Kerr metric in Cartesian Boyer-Linquist coordinate
void set_Kerr_metric_BL(const point_t &pos,
                        sym_tensor_rank2_spacetime &gKerrBL,
                        sym_tensor_rank2_spacetime (&d_gKerrBL)[4]){

constexpr double M_BH = 1.0; // Background BH metric mass
constexpr double r_sch = 2.*M_BH; // Schwarzschild radius in geometrized unit
double x = pos[0], y = pos[1], z = pos[2]; // short hand notation for spatial coordinates
double coords[4] = {0.0, x, y, z}; //General spacetime coordiantes
double r_floor = 1e-6; // Small floor value to avoid radial distance r goes to zero
double R = std::sqrt(x*x + y*y + z*z) + r_floor; // Radial distance that we will use TODO : Maybe not a great idea...
double R2 = R*R;
double R3 = R2*R;
double R4 = R2*R2;
double R5 = R2*R3;
// Define dimensionaless spin
// TODO : make it as parameter
const double J_ang = 0.1; //Angular momentum
const double a_ang = J_ang/M_BH; // Spin parameter in geometrical unit

// Some short hand notations
const double a2 = a_ang*a_ang;
double x2 = x*x, y2 = y*y, z2 = z*z;

// Some precomputations for Kerr

// Radial relation in terms of Cartesian component and spin
double r2 = (R2-a2 + std::sqrt((R2-a2)*(R2-a2) + 4.0*a2*z2))/2.0;
double r = std::sqrt(r2);

// Usual rho and Delta in Cartesian
double rho_sq = r2 + a2*z2/r2;
double DeltaKerr = r2-r_sch*r+a2;

// Azimuthal component of metric in Cartesian
double gphiphi = (r2+a2+r_sch*r*a2*(1.0-z2/r2)/rho_sq)*(1.0-z2/r2);
double gtphi = -r_sch*r*a_ang*(1.0-z2/r2)/rho_sq;

// tt component
gKerrBL(0,0) = -1.0+r_sch*r/rho_sq;

// tx (or xt) component
gKerrBL(0,1) = -y/(x2+y2)*gtphi;

// ty (or yt) component
gKerrBL(0,2) = x/(x2+y2)*gtphi;

// tz (or zt) component
gKerrBL(0,3) = 0.0;

// xx component
gKerrBL(1,1) = r2*x2/(rho_sq*DeltaKerr) + gphiphi*(y/(x2+y2))*(y/(x2+y2)) + x2*z2/(rho_sq*(r2-z2));

// yy component
gKerrBL(2,2) = r2*y2/(rho_sq*DeltaKerr) + gphiphi*(x/(x2+y2))*(x/(x2+y2)) + y2*z2/(rho_sq*(r2-z2));

// zz component
gKerrBL(3,3) = z2*(a2+r2)*(a2+r2)/(r2*rho_sq*DeltaKerr) 
             + rho_sq/(r2-z2)*(1.0-z2*(a2+r2)/(r2*rho_sq))*(1.0-z2*(a2+r2)/(r2*rho_sq));

// xy (or yx) component
gKerrBL(1,2) = r2*x*y/(rho_sq*DeltaKerr) - gphiphi*x*y/((x2+y2)*(x2+y2)) + x*y*z2/(rho_sq*(r2-z2));

// xz (or zx) component
gKerrBL(1,3) = (a2+r2)*x*z/(rho_sq*DeltaKerr) - x*z/(r2-z2)*(1.0-z2*(a2+r2)/(r2*rho_sq));

// yz (or zy) component
gKerrBL(2,3) = (a2+r2)*y*z/(rho_sq*DeltaKerr) - y*z/(r2-z2)*(1.0-z2*(a2+r2)/(r2*rho_sq));

// First derivative of metric

// tab components 
for(int i = 0; i < 4; ++i) {
  for(int j = 0; j < 4; ++j) {
    d_gKerrBL[0](i,j) = 0.0;
  }
}

//iab components

// Some more precomputations
double a2pr2 = a2+r2;
double delta1 = 1.0/DeltaKerr;
double sintheta_sq = 1.0-z2/r2;
double r21 = 1.0/r2;
double rho21 = 1.0/rho_sq;
double r2mz2 = r2-z2;
double r2mz21 = 1.0/r2mz2;
double x2py2 = x2+y2;
double x2py21= 1.0/x2py2;

double term1 = a2pr2*rho21;
double term2 = (2.0*a2pr2*a2pr2*rho21*rho21*z2*z)/(r2*r2);
 
// Derivative precomputations
double dr2dx = 2.0*r2*x/rho_sq;
double dr2dy = 2.0*r2*y/rho_sq;
double dr2dz = 2.0*a2pr2*z/rho_sq;

double drho2dx = 2.0*r21*rho21*x*(r2*r2-a2*z2);
double drho2dy = 2.0*r21*rho21*y*(r2*r2-a2*z2);
double drho2dz = 2.0*r21*r21*z*(a2*r2 + a2pr2*rho21*(r2*r2-a2*z2));

double ddeltadx = r*rho21*(2.0*r-r_sch)*x;
double ddeltady = r*rho21*(2.0*r-r_sch)*y;
double ddeltadz = (a2pr2*rho21*(2.0*r-r_sch)*z)/r;

double dgtphidx = (a_ang*rho21*rho21*r_sch*(drho2dx*r2mz2-x*(r2+z2)))/r;
double dgtphidy = (a_ang*rho21*rho21*r_sch*(drho2dy*r2mz2-y*(r2+z2)))/r;
double dgtphidz = a_ang*drho2dz*r*rho21*rho21*r_sch*sintheta_sq 
                - (a_ang*a2pr2*rho21*rho21*r_sch*sintheta_sq*z)/r + a_ang*r*rho21*r_sch*(2.0*r21*z - 2.0*a2pr2*r21*r21*rho21*z2*z);

double dsintheta_sqdx = 2.0*r21*rho21*x*z2;
double dsintheta_sqdy = 2.0*r21*rho21*y*z2;
double dsintheta_sqdz = -2.0*r21*z + 2.0*a2pr2*r21*rho21*rho21*z2*z;

double dgphiphidx = dsintheta_sqdx*r2 + dr2dx*sintheta_sq + a2*rho21*rho21*(dsintheta_sqdx*rho_sq*rho_sq 
                  - drho2dx*r*r_sch*sintheta_sq*sintheta_sq + r*r_sch*sintheta_sq*(2.0*dsintheta_sqdx*rho_sq 
                  + sintheta_sq*x));
double dgphiphidy = dsintheta_sqdy*r2 + dr2dy*sintheta_sq + a2*rho21*rho21*(dsintheta_sqdy*rho_sq*rho_sq 
                  - drho2dy*r*r_sch*sintheta_sq*sintheta_sq + r*r_sch*sintheta_sq*(2.0*dsintheta_sqdy*rho_sq 
                  + sintheta_sq*y));
double dgphiphidz = dsintheta_sqdz*r2 + dr2dz*sintheta_sq + (a2*a2*rho21*rho21*r_sch*sintheta_sq*sintheta_sq*z)/r 
                  + a2*rho21*rho21*(dsintheta_sqdz*rho_sq*rho_sq - drho2dz*r*r_sch*sintheta_sq*sintheta_sq 
                  + r*r_sch*sintheta_sq*(2.0*dsintheta_sqdz*rho_sq + sintheta_sq*z));

// Metric derivatives
// Derivatives with respect to x
d_gKerrBL[1](0,0) = r*rho21*rho21*r_sch*(x-drho2dx);
d_gKerrBL[1](0,1) = x2py21*(-dgtphidx + 2.0*gtphi*x*x2py21)*y;
d_gKerrBL[1](1,1) = -(ddeltadx*delta1*delta1*r2*rho21*x2) + delta1*r2*rho21*(2.0*x + 2.0*rho21*x2*x - drho2dx*rho21*x2) 
                  + x2py21*x2py21*(dgphiphidx - 4.0*gphiphi*x*x2py21)*y2 
                  - r2mz21*rho21*(-2.0*x + 2.0*r2*r2mz21*rho21*x2*x 
                  + drho2dx*rho21*x2)*z2;
d_gKerrBL[1](0,2) = x2py21*(gtphi + dgtphidx*x - 2.0*gtphi*x2*x2py21);
d_gKerrBL[1](1,2) = -(y*(ddeltadx*delta1*delta1*r2*rho21*x + delta1*r2*rho21*(-1.0 + drho2dx*rho21*x - 2.0*rho21*x2) 
                  + x2py21*x2py21*(gphiphi + dgphiphidx*x - 4.0*gphiphi*x2*x2py21) 
                  + r2mz21*rho21*(-1.0 + drho2dx*rho21*x 
                  + 2*r2*r2mz21*rho21*x2)*z2));
d_gKerrBL[1](2,2) = dgphiphidx*x2*x2py21*x2py21 + 2.0*gphiphi*x*x2py21*x2py21*(1.0 - 2.0*x2*x2py21) 
                   - rho21*y2*(delta1*r2*(ddeltadx*delta1 + rho21*(drho2dx - 2.0*x)) 
                   + r2mz21*rho21*(drho2dx + 2.0*r2*r2mz21*x)*z2);
d_gKerrBL[1](0,3) = 0.0;
d_gKerrBL[1](1,3) = -(z*(r2mz21 + delta1*term1*(-1.0 + ddeltadx*delta1*x + drho2dx*rho21*x) 
                  - 2.0*delta1*r2*rho21*rho21*x2 + r2mz21*(-2.0*rho21*rho21*x2 
                  + r21*term1*(-1.0 + rho21*x*(drho2dx + 2.0*x)))*z2 + 2.0*r2*r2mz21*r2mz21*rho21*x2
                  * (-1.0 + r21*term1*z2)));
d_gKerrBL[1](2,3) = -(y*z*(ddeltadx*delta1*delta1*term1 + delta1*rho21*(drho2dx*term1 - 2.0*r2*rho21*x) 
                  + r2mz21*rho21*(-2.0*r2*r2mz21*x - 2.0*rho21*x*z2 + r21*term1*(drho2dx + 2.0*(x + r2*r2mz21*x))*z2)));
d_gKerrBL[1](3,3) = -(a2pr2*a2pr2*delta1*drho2dx*r21*rho21*rho21*z2) - a2pr2*a2pr2*delta1*r21*rho21*(ddeltadx*delta1 
                  + 2.0*rho21*x)*z2 - 2.0*r2*r2mz21*r2mz21*x*(-1.0 + r21*term1*z2)*(-1.0 + r21*term1*z2) 
                  - drho2dx*r2mz21*(-1.0 + r21*term1*z2)*(1.0 + r21*(-1.0 + 2.0*rho_sq*rho21)*term1*z2) 
                  + 4.0*rho21*x*z2*(delta1*term1 + r2mz21*rho_sq*(rho21 - r21*term1)*(-1.0 + r21*term1*z2));

// Derivatives with respect to y
d_gKerrBL[2](0,0) = r*rho21*rho21*r_sch*(y-drho2dy);
d_gKerrBL[2](0,1) = x2py21*(-(dgtphidy*y) + gtphi*(-1.0 + 2.0*x2py21*y2));
d_gKerrBL[2](1,1) = -(ddeltady*delta1*delta1*r2*rho21*x2) - delta1*r2*rho21*rho21*x2*(drho2dy - 2.0*y) 
                  + x2py21*x2py21*(2.0*gphiphi*y - 4.0*gphiphi*x2py21*y2*y + dgphiphidy*y2) 
                  - r2mz21*rho21*rho21*x2*(drho2dy + 2.0*r2*r2mz21*y)*z2;
d_gKerrBL[2](0,2) = x*x2py21*(dgtphidy - 2.0*gtphi*x2py21*y);
d_gKerrBL[2](1,2) = -(x*(ddeltady*delta1*delta1*r2*rho21*y + delta1*r2*rho21*(-1.0 + drho2dy*rho21*y - 2.0*rho21*y2) 
                  + x2py21*x2py21*(gphiphi + dgphiphidy*y - 4.0*gphiphi*x2py21*y2) 
                  + r2mz21*rho21*(-1.0 + drho2dy*rho21*y + 2.0*r2*r2mz21*rho21*y2)*z2));
d_gKerrBL[2](2,2) = dgphiphidy*x2*x2py21*x2py21 + 2.0*delta1*r2*rho21*y - 4.0*gphiphi*x2*x2py21*x2py21*x2py21*y 
                  - delta1*r2*rho21*(-2.0*rho21*y2*y + ddeltady*delta1*y2 + drho2dy*rho21*y2) 
                  - r2mz21*rho21*(-2.0*y + 2.0*r2*r2mz21*rho21*y2*y + drho2dy*rho21*y2)*z2;
d_gKerrBL[2](0,3) = 0.0;
d_gKerrBL[2](1,3) = -(x*z*(ddeltady*delta1*delta1*term1 + delta1*rho21*(drho2dy*term1 - 2.0*r2*rho21*y) + r2mz21*rho21*(-2.0*r2*r2mz21*y 
                  - 2.0*rho21*y*z2 + r21*term1*(drho2dy + 2.0*(y + r2*r2mz21*y))*z2)));
d_gKerrBL[2](2,3) = -(z*(r2mz21 + delta1*term1*(-1.0 + ddeltady*delta1*y + drho2dy*rho21*y) 
                  - 2.0*delta1*r2*rho21*rho21*y2 
                  + r2mz21*(-2.0*rho21*rho21*y2 + r21*term1*(-1.0 + rho21*y*(drho2dy + 2.0*y)))*z2 
                  + 2.0*r2*r2mz21*r2mz21*rho21*y2
                  * (-1.0 + r21*term1*z2)));
d_gKerrBL[2](3,3) = -(a2pr2*a2pr2*delta1*drho2dy*r21*rho21*rho21*z2) - a2pr2*a2pr2*delta1*r21*rho21*(ddeltady*delta1 
                  + 2.0*rho21*y)*z2 
                  - 2.0*r2*r2mz21*r2mz21*y*(-1.0 + r21*term1*z2)*(-1.0 + r21*term1*z2) 
                  - drho2dy*r2mz21*(-1.0 + r21*term1*z2)*(1.0 + r21*(-1.0 + 2.0*rho_sq*rho21)*term1*z2) 
                  + 4.0*rho21*y*z2*(delta1*term1 + r2mz21*rho_sq*(rho21 - r21*term1)*(-1.0 + r21*term1*z2));

// Derivatives with respect to z
d_gKerrBL[3](0,0) = -drho2dz*r*rho21*rho21*r_sch + (rho21*r_sch*term1*z)/r;
d_gKerrBL[3](0,1) = -dgtphidz*x2py21*y;
d_gKerrBL[3](1,1) = -ddeltadz*delta1*delta1*r2*rho21*x2 - delta1*drho2dz*r2*rho21*rho21*x2 + dgphiphidz*x2py21*x2py21*y2 
                  + 2.0*r2mz21*rho21*x2*z + 2.0*delta1*rho21*term1*x2*z - drho2dz*r2mz21*rho21*rho21*x2*z2 
                  - r2mz21*r2mz21*rho21*x2*(-2.0*z + 2.0*term1*z)*z2;
d_gKerrBL[3](0,2) = dgtphidz*x*x2py21;
d_gKerrBL[3](1,2) = -ddeltadz*delta1*delta1*r2*rho21*x*y - delta1*drho2dz*r2*rho21*rho21*x*y - dgphiphidz*x*x2py21*x2py21*y 
                  + 2.0*r2mz21*rho21*x*y*z + 2.0*delta1*rho21*term1*x*y*z - drho2dz*r2mz21*rho21*rho21*x*y*z2 
                  - r2mz21*r2mz21*rho21*x*y*(-2.0*z + 2.0*term1*z)*z2;
d_gKerrBL[3](2,2) = dgphiphidz*x2*x2py21*x2py21 - ddeltadz*delta1*delta1*r2*rho21*y2 - delta1*drho2dz*r2*rho21*rho21*y2 
                  + 2.0*r2mz21*rho21*y2*z + 2.0*delta1*rho21*term1*y2*z - drho2dz*r2mz21*rho21*rho21*y2*z2 
                  - r2mz21*r2mz21*rho21*y2*(-2.0*z + 2.0*term1*z)*z2;
d_gKerrBL[3](0,3) = 0.0;
d_gKerrBL[3](1,3) = delta1*term1*x - ddeltadz*delta1*delta1*term1*x*z 
                  - delta1*drho2dz*rho21*term1*x*z + 2*delta1*rho21*term1*x*z2 
                  - r2mz21*x*(1.0 - r21*term1*z2) + r2mz21*r2mz21*x*z*(-2.0*z + 2.0*term1*z)*(1.0 - r21*term1*z2) 
                  - r2mz21*x*z*(term2 - 2.0*r21*term1*z - 2.0*r21*rho21*term1*z2*z + drho2dz*r21*rho21*term1*z2);
d_gKerrBL[3](2,3) = delta1*term1*y - ddeltadz*delta1*delta1*term1*y*z - delta1*drho2dz*rho21*term1*y*z 
                  + 2*delta1*rho21*term1*y*z2 
                  - r2mz21*y*(1.0 - r21*term1*z2) + r2mz21*r2mz21*y*z*(-2.0*z + 2.0*term1*z)*(1.0 - r21*term1*z2) 
                  - r2mz21*y*z*(term2 - 2.0*r21*term1*z - 2.0*r21*rho21*term1*z2*z + drho2dz*r21*rho21*term1*z2);
d_gKerrBL[3](3,3) = 2.0*a2pr2*a2pr2*delta1*r21*rho21*z - (2.0*a2pr2*a2pr2*a2pr2*delta1*rho21*rho21*z2*z)/(r2*r2) 
                  + 4.0*a2pr2*a2pr2*delta1*r21*rho21*rho21*z2*z - a2pr2*a2pr2*ddeltadz*delta1*delta1*r21*rho21*z2 
                  - a2pr2*a2pr2*delta1*drho2dz*r21*rho21*rho21*z2 
                  + drho2dz*r2mz21*(1.0 - r21*term1*z2)*(1.0 - r21*term1*z2) 
                  - r2mz21*r2mz21*rho_sq*(-2.0*z + 2.0*term1*z)*(1.0 - r21*term1*z2)*(1.0 - r21*term1*z2) 
                  + 2.0*r2mz21*rho_sq*(1.0 - r21*term1*z2)*(term2 - 2.0*r21*term1*z - 2.0*r21*rho21*term1*z2*z 
                  + drho2dz*r21*rho21*term1*z2);
} // Kerr in BL


//Metric from RNSID
void set_RNSID_metri(){

  //TODO : Finsh this later or might not be useful
  log_one(info)<<"Sorry this is not implemented yet"<<std::endl;
  assert(false);
  //#include "rnsid.h"
} // RNSID

} //namespace backgroud_metric


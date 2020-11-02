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
//#include "params.h"
#include "density_profiles.h"
void set_TOV_metric(const point_t & pos,
                    sym_tensor_rank2 & gTOV,
                    sym_tensor_rank2 & inv_gTOV,
                    sym_tensor_rank2 (&d_gTOV)[4]){
  double alpha2, dalpha2dr, beta2, dbeta2dr;
  double x = pos[0], y = pos[1], z = pos[2]; // short hand notation for spatial coordinates
  double r = std::sqrt(x*x + y*y + z*z) + 1.0e-7; // Radial distance that we will use TODO : singular
  double coords[4] = {0.0, x, y, z};
  double r2 = r*r;
  double r3 = r*r2;
  double r4 = r2*r2;
  double r5 = r2*r3;
  //sym_tensor_rank2 dt_gTOV{0};
  //sym_tensor_rank2 dx_gTOV{0};
  //sym_tensor_rank2 dy_gTOV{0};
  //sym_tensor_rank2 dz_gTOV{0};
  
  for(int mu = 0; mu < 4; mu++){
    for(int nu = 0; nu < 4; nu++){
      (d_gTOV[0])(mu,nu) = 0;
    }
  }

  // In the interior of the star: TOV metric
  if(r < param::sphere_radius){
    alpha2 = density_profiles::spherical_alpha2(r/param::sphere_radius);
    dalpha2dr = density_profiles::spherical_dalpha2_dr(r/param::sphere_radius);
    beta2 = density_profiles::spherical_beta2(r/param::sphere_radius);
    dbeta2dr = density_profiles::spherical_dbeta2_dr(r/param::sphere_radius);
    //log_one(info) << "alpha: "<<alpha2 <<"beta:" << beta2<< std::endl;
  } else {
    double r_sch = 1 - 2*GNEWT*param::sphere_mass/(C_LIGHT_CGS*C_LIGHT_CGS);
    alpha2 = (1 - r_sch / r);
    dalpha2dr = 2*r_sch/r2;
    beta2 = r/ (r - r_sch);
    dbeta2dr = -r_sch / (r - r_sch);
    //log_one(info) << "alpha: "<<alpha2 <<"beta:" << beta2<< std::endl;
  }
  gTOV(0,0) = -alpha2;
  inv_gTOV(0,0) = -1.0/alpha2;
  //dx_gTOV(0,0) = -dalpha2dr * coords[1] / r; 
  //dy_gTOV(0,0) = -dalpha2dr * coords[2] / r;
  //dz_gTOV(0,0) = -dalpha2dr * coords[3] / r;
  for(int i = 1; i< 4; i++) {
    (d_gTOV[i])(0,0) = -dalpha2dr * coords[i] / r;
  }
  //log_one(info) << "xyz: " << coords[1] << coords[2] << coords[3] << std::endl;
  for(int i = 1; i < 4; i++){
    for(int j = 1; j < 4; j++){
      gTOV(i,j) = (int)(i==j) + coords[i] * coords[j] * (beta2-1) / r2;
      inv_gTOV(i,j) = (int)(i==j) + coords[i] * coords[j] * (1.0/beta2-1) / r2;
      //log_one(info) << "g"<< i <<j << " : "<< (int)(i==j) << gTOV(i,j)<<std::endl;
      //product rule
      /*
      dx_gTOV(i,j) = ((i==1)*coords[j]+(j==1)*coords[i]) / r2 * (beta2 - 1)
                    + coords[i]*coords[j] * (-2*coords[1]/r4) * (beta2 - 1)
                    + coords[i]*coords[j] / r2 * (dbeta2dr);
      dy_gTOV(i,j) = ((i==2)*coords[j]+(j==2)*coords[i]) / r2 * (beta2 - 1)
                    + coords[i]*coords[j] * (-2*coords[2]/r4) * (beta2 - 1)
                    + coords[i]*coords[j] / r2 * (dbeta2dr);
      dz_gTOV(i,j) = ((i==3)*coords[j]+(j==3)*coords[i]) / r2 * (beta2 - 1)
                    + coords[i]*coords[j] * (-2*coords[3]/r4) * (beta2 - 1)
                    + coords[i]*coords[j] / r2 * (dbeta2dr);
      */
      for(int k = 1; k < 4; ++k) {
        (d_gTOV[k])(i,j) = ((int)(i==k) * coords[j] + (int)(j==k)*coords[i]) / r2 * (beta2 - 1)
                           + coords[i]*coords[j] * (-2*coords[k]/r4) * (beta2 - 1)
                           + coords[i]*coords[j] / r2 * (dbeta2dr);
      }
    }
  }
}

// HL : Adding/fixing now
#if 0
// Flat Minkowski metric in Carteisan coordinate
void set_Minkowski_metric(sym_tensor_rank2 & gMinkowski,
                          sym_tensor_rank2 (&d_gMinkowski)[4]){

 // Rank = 2, Dim = 4 -> 10 independent quantities
 gMinkowski(0,0) = -1.0; //tt
 gMinkowski(0,1) =  0.0; //tx
 gMinkowski(0,2) =  0.0; //ty
 gMinkowski(0,3) =  0.0; //tz
 gMinkowski(1,1) =  1.0; //xx
 gMinkowski(1,2) =  0.0; //xy
 gMinkowski(1,3) =  0.0; //xz
 gMinkowski(2,2) =  1.0; //yy
 gMinkowski(2,3) =  0.0; //yz
 gMinkowski(3,3) =  1.0; //zz

 // First derivative of metric
  for(int i = 0; i < 4; ++i) {
    for(int j = 0; j < 4; ++j) {
      for(int k = 0; k < 4; ++j){
        d_gMinkowski[i](j,k) = 0.0;
    }
   }
  }
} // Minkowski

//TODO : Check BH backgroud metics into Sparta to see Ham and mom constraint values
//TODO : Put metrics in Boyer-Lindquist coordinate form too. 
// Static spherically symmetric metric (i.e. Schwarzschild) in Kerr-Schild coordinate
void set_Schwarzschild_metric_KS(const point_t &pos,
                                 sym_tensor_rank2 & gSchwarzKS,
                                 sym_tensor_rank2 (&d_gSchwarzKS)[4]){

const double M_back = 1.0;
//const double r_sch = 2.*M_back*gc/(C_LIGHT_CGS*C_LIGHT_CGS); // Schwarzschild radiusa in cgs
const double r_sch = 2.*M_back; // Schwarzschild radius in geometrical unit
double coords[4] = {0}; //General spacetime coordiantes
double t = coords[0]; // short hand notation for time coordinate
double x = coords[1], y = coords[2], z = coords[3]; // short hand notation for spatial coordinates
double r_real = std::sqrt(x*x + y*y + z*z); // this is true radial distance
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

gSchwarz(0,0) = -alpha*alpha + beta_sum; //tt

//tx, ty, tz components
for(int i = 1; i < 4; ++i) {
  gSchwarz(0,i) = beta_d[i];
}

// ij (spatial) components
for(int i = 1; i < 4; ++i) {
  for(int j = i; j < 4; ++j) {
    gSchwarz(i,j) = (i==j) + r_sch*coords[i]*coords[j]/(r3);
  }
}

//First derivative of metric. 

sym_tensor_rank2 dt_gSchwarz{0}; // partial_t g_ab
sym_tensor_rank2 dx_gSchwarz{0}; // partial_x g_ab
sym_tensor_rank2 dy_gSchwarz{0}; // partial_y g_ab
sym_tensor_rank2 dz_gSchwarz{0}; // partial_z g_ab

//tab components. All zeros
for(int i = 0; i < 4; ++i) {
  for(int j = i; j < 4; ++j) {
    dt_gSchwarz(i,j) = 0.0;
  }
}

//itj (or ijt) components
for(int j = 1; j < 4; ++j){
    dx_gSchwarz(0,j) = r_sch*(1==j)/(r*(r+r_sch));
    dy_gSchwarz(0,j) = r_sch*(2==j)/(r*(r+r_sch));
    dz_gSchwarz(0,j) = r_sch*(3==j)/(r*(r+r_sch));
}

//ijk (all spatial) components
  for(int j = 1; j < 4; ++j) {
    for(int k = j; k < 4; ++k) {
      dx_gSchwarz(j,k) = r_sch/(r3)*((1==k) + (j==k)) 
                           + 3*r_sch*coords[1]*coords[j]*coords[k]/r5;
      dy_gSchwarz(j,k) = r_sch/(r3)*((2==k) + (j==k)) 
                           + 3*r_sch*coords[2]*coords[j]*coords[k]/r5;
      dz_gSchwarz(j,k) = r_sch/(r3)*((3==k) + (j==k)) 
                           + 3*r_sch*coords[3]*coords[j]*coords[k]/r5;
    }
  }

} // Schwarzschild in KS

// Schwarzschild in isotropic coordinates
void set_Schwarzschild_metric_BL(const point_t &pos,
                                 sym_tensor_rank2 & gSchwarzIso,
                                 sym_tensor_rank2 (&d_gSchwarzIso)[4]){

} // Schwarzschild in BL

// Axisymmetic metric (i.e. Kerr) in Kerr-Schild Cartesian coordinate
// NOTE : I keep both Schwarzschild and Kerr for now for sanity check. 
//        Once everything looks fine, I will remove Schwarzschild since
//        a->0 (or J->0) in Kerr will return Schwarzschild
void set_Kerr_metric_KS(const point_t &pos,
                        sym_tensor_rank2 & gKerrKS,
                        sym_tensor_rank2 (&d_gKerrKS)[4]){

sym_tensor_rank2 gKerr{0};

//Define dimensionaless spin
//TODO : make it as parameter
const double J_ang = 0.1; //Angular momentum
//const double a_ang = J_ang/(M_back*C_LIGHT_CGS); // Spin parameter in cgs
const double a_ang = J_ang/M_back; // Spin parameter in geometrical unit

//Some short hand notation
const double a2 = a_ang*a_ang;
double x2 = x*x, y2 = y*y, z2 = z*z;

//Define scalar quantities 
//double f_scalar = 2.*gc*M_back*r3/(r4+a2*z2);
double f_scalar = 2.*gc*r3/(r4+a2*z2); // in geometrical unit

//Define k 4-vector in covariant form
double k_vec[4];
k_vec[0] = 1.0;
k_vec[1] = (r*x+a_ang*y)/(r2+a2); 
k_vec[2] = (r*y-a_ang*x)/(r2+a2); 
k_vec[3] = z/r; 

for(int i = 0; i < 4; ++i){
  for(int j = 0; j < 4; ++j){
    gKerr(i,j) = gMinkowski(i,j) + f_scalar*k_vec[i]*k_vec[j];
  }
}

//First derivative of metric
sym_tensor_rank2 dt_gKerr{0}; // partial_t g_ab
sym_tensor_rank2 dx_gKerr{0}; // partial_x g_ab
sym_tensor_rank2 dy_gKerr{0}; // partial_y g_ab
sym_tensor_rank2 dz_gKerr{0}; // partial_z g_ab


//Define derivative qunatities
//Derivatives of scalar
double d_f[4];
d_f[0] = 0.0;
#if 0 //HL : keep this until we have correct unit conversion
d_f[1] = 6.0*gc*M_back*x*r/(a2*z2+r4)-8.0*gc*M_back*x*r5/((a2*z2+r4)*(a2*z2+r4));
d_f[2] = 6.0*gc*M_back*y*r/(a2*z2+r4)-8.0*gc*M_back*y*r5/((a2*z2+r4)*(a2*z2+r4));
d_f[3] = 6.0*gc*M_back*y*r/(a2*z2+r4)-2.0*gc*M_back*r3*(2.0*a2*z+4.0*z*r2)/((a2*z2+r4)*(a2*z2+r4));
#endif
d_f[1] = 6.0*M_back*x*r/(a2*z2+r4)-8.0*M_back*x*r5/((a2*z2+r4)*(a2*z2+r4));
d_f[2] = 6.0*M_back*y*r/(a2*z2+r4)-8.0*M_back*y*r5/((a2*z2+r4)*(a2*z2+r4));
d_f[3] = 6.0*M_back*y*r/(a2*z2+r4)-2.0*M_back*r3*(2.0*a2*z+4.0*z*r2)/((a2*z2+r4)*(a2*z2+r4));
//Derivatives of k-vector
double d_k[4][4];
//ta or at components are zero
for(int i = 1; i < 4; ++i){
   d_k[0][0] = 0.0;
   d_k[0][i] = 0.0;
   d_k[i][0] = 0.0;

}
//ij (spatial) components
d_k[1][1] = (x2/r+r)/(a2+r2)-2*x*(a_ang*y+x*r)/((a2+r2)*(a2+r2));
d_k[1][2] = (x*y/r+a_ang)/(a2+r2)-2*y*(a_ang*y+x*r)/((a2+r2)*(a2+r2));
d_k[1][3] = x*z/(r*(a2+r2))-2*z*(a_ang*y+x*r)/((a2+r2)*(a2+r2));
d_k[2][1] = (x*y/r-a_ang)/(a2+r2)-2*x*(y*r-a_ang*x)/((a2+r2)*(a2+r2));
d_k[2][2] = (y2/r+r)/(a2+r2)-2*y*(y*r-a_ang*x)/((a2+r2)*(a2+r2));
d_k[2][3] = y*z/(r*(a2+r2))-2*z*(y*r-a_ang*x)/((a2+r2)*(a2+r2));
d_k[3][1] = -x*z/r3;
d_k[3][2] = -y*z/r3;
d_k[3][3] = -z2/r3 + 1/r;

//Define first derivative of Kerr metric
  for(int j = 0; j < 4; ++j) {
    for(int k = j; k < 4; ++k) {
      dt_gKerr(j,k) = d_f[0]*k_vec[j]*k_vec[k] 
                      + f_scalar*(d_k[k][0]*k_vec[j] + k_vec[0]*d_k[k][j]);
      dx_gKerr(j,k) = d_f[1]*k_vec[j]*k_vec[k] 
                      + f_scalar*(d_k[k][1]*k_vec[j] + k_vec[1]*d_k[k][j]);
      dy_gKerr(j,k) = d_f[2]*k_vec[j]*k_vec[k] 
                      + f_scalar*(d_k[k][2]*k_vec[j] + k_vec[2]*d_k[k][j]);
      dz_gKerr(j,k) = d_f[3]*k_vec[j]*k_vec[k] 
                      + f_scalar*(d_k[k][3]*k_vec[j] + k_vec[3]*d_k[k][j]);
    }
  }
} // Kerr in KS

// Kerr metric in Boyer-Linquist coordinate
void set_Kerr_metric_BL(const point_t &pos,
                        sym_tensor_rank2 &gKerrBL,
                        sym_tensor_rank2 (&d_gKerrBL)[4]){

} // Kerr in BL

//Metric from RNSID
void set_RNSID_metri(){

  //TODO : Finsh this later or might not be useful

  //#include "rnsid.h"
} // RNSID
#endif

} //namespace backgroud_metric


//#endif

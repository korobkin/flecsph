/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2018 Triad National Security, LLC
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
 * @file fmm.h
 * @brief Functions used in the FMM computation
 */

#pragma once

#include "params.h"
#include "tree.h"

namespace fmm {
  using namespace param;

  int index(int i, int j, int k, int l){
    return ((i*3+j)*3+k)*3+l;
  }

  int index(int i, int j, int k){
    return (i*3+j)*3+k;
  }

  int index(int i, int j){
    return i*3+j;
  }

  void compute_momentum(
    std::array<double,91>& X,
    std::array<double,27>& H, 
    std::array<double,9>& Q,
    std::vector<body*> bs, 
    const point_t& cofm_center)
  {
    for(int b = 0 ; b < bs.size(); ++b){  
      double mb = bs[b]->mass();
      point_t q = bs[b]->coordinates()-cofm_center;
      double q2 = q[0]*q[0]+q[1]*q[1]+q[2]*q[2];
      double q4 = q2*q2; 
      for(int i = 0 ; i < 3; ++i){
        for(int j = 0 ; j < 3; ++j){
          Q[index(i,j)] += mb*(3*q[i]*q[j]-(i==j)*q2);
          for(int k = 0 ; k < 3; ++k){
            H[index(i,j,k)] += mb*(
              15*q[i]*q[j]*q[k]
              -3*q2*((i==j)*q[k]+(j==k)*q[i]+(i==k)*q[j])
            );
            for(int l = 0 ; l < 3; ++l){
              X[index(i,j,k,l)] += mb*(
                105*q[i]*q[j]*q[k]*q[l]-
                  15*q2*(
                    (i==j)*q[k]*q[l]+(i==l)*q[j]*q[k]+(i==k)*q[j]*q[l]+
                    (j==l)*q[i]*q[k]+(j==k)*q[i]*q[l]+(l==k)*q[i]*q[j]
                  )+
                  3*q4*((i==j)*(k==l)+(i==k)*(j==l)+(i==l)*(j==k))
              );
            }
          }
        }
      }
    }
  }

  /*
  * @brief Compute the gravitation interation
  * Return the computed value if needed for direct particle interaction
  *
  * The Sink is the one on which I compute the fc, dfcdr, dfcdrdr
  * The Source is the distant box
  */
  inline
  point_t gravitation_p2p(
    point_t & fc,
    const point_t& sink_coordinates,
    const point_t& source_coordinates,
    const double& sm)
  {
    double dist = flecsi::distance(sink_coordinates,source_coordinates);
    point_t res = -gravitational_constant*sm/(dist*dist*dist)*
      (sink_coordinates-source_coordinates);
    fc += res;
    return res;
  }


  /*
  * @brief Compute the gravitation interation
  * Return the computed value if needed for direct particle interaction
  *
  * The Sink is the one on which I compute the fc, dfcdr, dfcdrdr
  * The Source is the distant box
  */
  inline
  void gravitation_fc(
    point_t & fc,
    const point_t& local_coordinates,
    const point_t& dist_coordinates,
    const double& sm, 
    const std::array<double,91>& X, 
    const std::array<double,27>& H, 
    const std::array<double,9>& Q)
  {
    double d = flecsi::distance(local_coordinates,dist_coordinates);
    double d2 = d*d; 
    double d3 = d2*d; 
    double d5 = d3*d2; 
    double d7 = d5*d2;
    double d9 = d7* d2; 
    double d11 = d9*d2;
    point_t r = dist_coordinates-local_coordinates; 

    for(int m = 0; m < 3; ++m){
      // Monopole
      fc -= gravitational_constant*sm*r[m]/d3;
      for(int i = 0 ; i < 3; ++i){
        // Quadrupole 
        fc += Q[index(i,m)]*r[i]/d5;
        for(int j = 0 ; j < 3; ++j){
          // Quadrupole 
          fc -= 5./2.*Q[index(i,j)]*r[i]*r[j]*r[m]/d7;
          // Octopole 
          fc += .5*H[index(i,j,m)]*r[i]*r[j]/d7; 
          for(int k = 0 ; k < 3; ++k){
            // Octopole 
            fc -= 7./6.*H[index(i,j,k)]*r[i]*r[j]*r[k]*r[m]/d9;
            // Hexadecapole 
            fc += 1./6.*X[index(i,j,k,m)]*r[i]*r[j]*r[k]/d9;  
            for(int l = 0; l < 3; ++l){
              // Hexadecapole 
              fc -= 9./24.*X[index(i,j,k,l)]*r[i]*r[j]*r[k]*r[l]*r[m]/d11;
            }
          }
        }
      }
    }
    #if 0 
    double dist = flecsi::distance(sink_coordinates,source_coordinates);
    point_t res = -gravitational_constant*source_mass/(dist*dist*dist)*
      (sink_coordinates-source_coordinates);
    fc += res;
    return res;
    #endif 
  }

  /*
  * @brief Compute the Jacobian (dfcdr) matrix on the Sink from the Source
  */
  inline
  void gravitation_dfcdr(
    double res[9],
    const point_t& local_coordinates,
    const point_t& dist_coordinates,
    const double& sm, 
    const std::array<double,91>& X, 
    const std::array<double,27>& H, 
    const std::array<double,9>& Q)
  {
    double d = flecsi::distance(local_coordinates,dist_coordinates);
    double d2 = d*d; 
    double d3 = d2*d; 
    double d5 = d3*d2; 
    double d7 = d5*d2;
    double d9 = d7* d2; 
    double d11 = d9*d2;
    point_t r = dist_coordinates-local_coordinates; 

    for(int m = 0; m < 3; ++m){
      for(int q = 0; q < 3; ++q){
        int pr = m*3+q;
        // Monopole
        res[pr] += gravitational_constant*sm/d3*(3*r[m]*r[q]/d2-(q==m));
        // Quadrupole
        res[pr] += Q[index(m,q)]/d5; 
        for(int i = 0 ; i < 3; ++i){
          // Quadrupole
          res[pr] -= 5*(Q[index(i,m)]*r[q]+Q[index(i,q)]*r[m])*r[i]/d7; 
          // Octopole 
          res[pr] += H[index(i,q,m)]*r[i]/d7;
          for(int j = 0 ; j < 3; ++j){
            // Quadrupole 
            res[pr] += (35./2.*r[m]*r[q]/d2-5./2.*(m==q))*
              Q[index(i,j)]*r[i]*r[j]/d7; 
            // Octopole 
            res[pr] -= 7./2.*(H[index(i,j,m)]*r[q]+
              H[index(i,j,q)]*r[m])*r[i]*r[j]/d9;
            // Hexadecapole
            res[pr] += .5*X[index(i,j,q,m)]*r[i]*r[j]/d9; 
            for(int k = 0 ; k < 3; ++k){
              // Octopole 
              res[pr] += 7./6.*(9.*r[m]*r[q]/d2-(q==m))*
                H[index(i,j,k)]*r[i]*r[j]*r[k]/d9; 
              // Hexadecapole 
              res[pr] -= 9./6.*(X[index(i,j,k,m)]*r[q]+
                X[index(i,j,k,q)]*r[m])*r[i]*r[j]*r[k]/d11;
              for(int l = 0; l < 3; ++l){
                // Hexadecapole 
                res[pr] += 9./24.*(11.*r[m]*r[q]/d2-(q==m))*
                  X[index(i,j,k,l)]*r[i]*r[j]*r[k]*r[l]/d11;
              }
            }
          }
        }
      }
    }
#if 0 
    double dist = flecsi::distance(sink_coordinates,source_coordinates);
    double dist_2 = dist*dist;
    point_t diffPos =  sink_coordinates - source_coordinates;
    double jacobicoeff = -gravitational_constant*source_mass/(dist_2*dist);
    for(int i = 0; i < 9; ++i){
      int a = i/3; int b = i%3;
      double valjacobian = jacobicoeff*((a==b)-3*diffPos[a]*diffPos[b]/(dist_2));
      dfcdr[i] += valjacobian;
    }
#endif 

  }

  /*
  * @brief Compute the Hessian (dfcdrdr) matrix on the Sink from the Source
  */
  inline
  void gravitation_dfcdrdr(
    double res[27],
    const point_t& local_coordinates,
    const point_t& dist_coordinates,
    const double& sm, 
    const std::array<double,91>& X, 
    const std::array<double,27>& H, 
    const std::array<double,9>& Q)
  {

    double d = flecsi::distance(local_coordinates,dist_coordinates);
    const double d2 = d*d; 
    const double d3 = d2*d; 
    const double d5 = d3*d2; 
    const double d7 = d5*d2;
    const double d9 = d7*d2; 
    const double d11 = d9*d2;
    const double d13 = d11*d2; 
    const double d15 = d13*d2; 
    point_t r = dist_coordinates-local_coordinates; 

    for(int m = 0; m < 3; ++m){
      for(int q = 0; q < 3; ++q){
        for(int s = 0 ; s < 3; ++s){
          int pr = (m*3+q)*3+s;
          // Monopole 
          res[pr] += gravitational_constant*3*sm/d5*(
            (m==q)*r[s]+(q==s)*r[m]+(m==s)*r[q]-5*r[m]*r[q]*r[s]/d2);
          // Quadrupole
          res[pr] += 5./d7*
            (Q[index(m,q)]*r[s]+Q[index(s,m)]*r[q]+Q[index(s,q)]*r[m]);
          // Octopole 
          res[pr] += H[index(s,q,m)]/d7;
          for(int i = 0 ; i < 3; ++i){
            // Quadrupole 
            res[pr] -= Q[index(i,m)]*(q==s)+
                       Q[index(i,q)]*(m==s)+
                       Q[index(i,s)]*(m==q)*5.*r[i]/d7;
            res[pr] += Q[index(i,m)]*r[q]*r[s]+Q[index(i,q)]*r[m]*r[s]+
              Q[index(i,s)]*r[m]*r[q]*35.*r[i]/d9; 
            // Octopole 
            res[pr] -= (H[index(i,q,m)]*r[s]+
                        H[index(i,s,m)]*r[q]+
                        H[index(i,s,q)]*r[m])*7.*r[i]/d9;
            // Hexadecapole 
            res[pr] += X[index(i,s,q,m)]*r[i]/d9; 
            for(int j = 0 ; j < 3; ++j){
              // Quadrupole
              res[pr] += 35./2.*((m==q)*r[s]+
                                 (m==s)*r[q]+
                                 (s==q)*r[m])*Q[index(i,j)]*r[i]*r[j]/d9;
              res[pr] += 315./2.*Q[index(i,j)]*r[i]*r[j]*r[m]*r[q]*r[s]/d11;
              // Octopole 
              res[pr] -= 7./2.*(H[index(i,j,m)]*(q==s)+
                                H[index(i,j,q)]*(m==s)+
                                H[index(i,j,s)]*(q==m))*r[i]*r[j]/d9;
              res[pr] += 63./2.*(H[index(i,j,m)]*r[q]*r[s]+
                                 H[index(i,j,q)]*r[m]*r[s]+
                                 H[index(i,j,s)]*r[m]*r[q])*r[i]*r[j]/d11;
              // Hexadecapole 
              res[pr] -= 9./2.*(X[index(i,j,q,m)]*r[s]+
                                X[index(i,j,s,m)]*r[q]+
                                X[index(i,j,s,q)]*r[m])*r[i]*r[j]/d11;
              for(int k = 0 ; k < 3; ++k){
                // Octopole 
                res[pr] += 63./6.*((q==m)*r[s]+
                                   (m==s)*r[q]+
                                   (q==s)*r[m])*
                                   H[index(i,j,k)]*r[i]*r[j]*r[k]/d11;
                res[pr] -= 693./6.*H[index(i,j,k)]*
                  r[i]*r[j]*r[k]*r[m]*r[q]*r[s]/d13;
                // Hexadecapole 
                res[pr] -= 9./6.*(X[index(i,j,k,m)]*(q==s)+
                                  X[index(i,j,k,q)]*(s==m)+
                                  X[index(i,j,k,s)]*(q==m))*r[i]*r[j]*r[k]/d11;
                res[pr] += 99./6.*(X[index(i,j,k,m)]*r[q]*r[s]+
                                   X[index(i,j,k,q)]*r[m]*r[s]+
                                   X[index(i,j,k,s)]*r[m]*r[q])*r[i]*r[j]*r[k]/d13;
                for(int l = 0; l < 3; ++l){
                  // Hexadecapole
                  res[pr] += 99./24.*((q==m)*r[s]+
                                      (m==s)*r[q]+
                                      (q==s)*r[m])*
                                      X[index(i,j,k,l)]*r[i]*r[j]*r[k]*r[l]/d13;
                  res[pr] -= 1287./24.*X[index(i,j,k,l)]*
                      r[i]*r[j]*r[k]*r[l]*r[m]*r[q]*r[s]/d15; 
                }
              }
            }
          }
        }
      }
    }


#if 0
    double dist = flecsi::distance(sink_coordinates,source_coordinates);
    double dist_2 = dist*dist;
    point_t diffPos =  sink_coordinates - source_coordinates;
    double hessiancoeff = -gravitational_constant*3.0*source_mass/(dist_2*dist_2*dist);
    for(int i = 0 ; i < 27 ; ++i){
      int a = i/9; int b = (i%9)/3; int c = i%3;
      double term_1 = (a==b)*diffPos[c]+(c==a)*diffPos[b]+(b==c)*diffPos[a];
      double valhessian = hessiancoeff *
        ( 5.0/(dist_2)*diffPos[a]*diffPos[b]*diffPos[c] - term_1) ;
      dfcdrdr[i] += valhessian;
    }
#endif 
  }

  /*
  * @brief Taylor expansion of degree 2 using the computed function, jacobi,
  * hessian and the targeted particle
  */
  void interation_c2p(
    point_t& fc,
    double dfcdr[9],
    double dfcdrdr[27],
    point_t cofm_coordinates,
    body* sink)
  {
    point_t part_coordinates = sink->coordinates();

    point_t diffPos = part_coordinates - cofm_coordinates;
    point_t grav = fc;
    // The Jacobi
    for(size_t i=0;i<gdimension;++i){
      for(size_t j=0;j<gdimension;++j){
        grav[i] += dfcdr[i*gdimension+j]*diffPos[j];
      } // for
    } // for
    // The hessian
    double tmpMatrix[gdimension*gdimension] = {};
    for(size_t i=0;i<gdimension;++i){
      for(size_t j=0;j<gdimension;++j){
        for(size_t k=0;k<gdimension;++k){
          tmpMatrix[i*gdimension+j] +=
            diffPos[k]*dfcdrdr[i*gdimension*gdimension+j*gdimension+k];
        } // for
      } // for
    } // for
    double tmpVector[gdimension] = {};
    for(size_t i=0;i<gdimension;++i){
      for(size_t j=0;j<gdimension;++j){
        tmpVector[j] += tmpMatrix[i*gdimension+j]*diffPos[i];
      } // for
    } // for
    for(size_t i=0;i<gdimension;++i){
      grav[i] += 0.5*tmpVector[i];
    } // for
    sink->setAcceleration(grav+sink->getAcceleration());
  }

} // namespace fmm


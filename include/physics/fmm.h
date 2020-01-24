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

#define QUAD
#define OCTO
//#define HEXA

#pragma once

#include "params.h"
#include "tree.h"

namespace fmm {
  using namespace param;


  struct fmm_comms{
    fmm_comms(){}
    fmm_comms(
      const point_t& _coords, 
      const double& _T,
      const tensor_u<double, symmetry_type::symmetric, 3, 3, 3, 3>& _X,
      const tensor_u<double, symmetry_type::symmetric, 3, 3, 3>& _H, 
      const tensor_u<double, symmetry_type::symmetric, 3, 3>& _Q)
    : coords(_coords), T(_T), X(_X), H(_H), Q(_Q){}
    point_t coords; 
    double T; 
    tensor_u<double, symmetry_type::symmetric, 3, 3, 3, 3> X; 
    tensor_u<double, symmetry_type::symmetric, 3, 3, 3> H; 
    tensor_u<double, symmetry_type::symmetric, 3, 3> Q; 
  };


  /**
  * Compute momenta for a specific q (vector between the two
  *  entities, mass of entity 0 and mass of entity 1 using 
  * the reduced mass. 
  */
  inline double compute_XHQ(
    tensor_u<double, symmetry_type::symmetric, 3, 3, 3, 3>& X,
    tensor_u<double, symmetry_type::symmetric, 3, 3, 3>& H,
    tensor_u<double, symmetry_type::symmetric, 3, 3>& Q,
    const double& m0, const double& m1, const point_t& q)
  {
    const double m = m0*m1/(m0+m1);
    const double q2 = q[0]*q[0]+q[1]*q[1]+q[2]*q[2]; 
    const double q4 = q2*q2; 
    for(int i = 0 ; i < 3; ++i){
      for(int j = i ; j < 3; ++j){
        // Quadrupole 
        Q(i,j) += m*(3.*q[i]*q[j]-(i==j)*q2); 
        for(int k = j ; k < 3; ++k){
          // Octopole
          H(i,j,k) += m*(
            15.*q[i]*q[j]*q[k]
            -3.*q2*((i==j)*q[k]+(j==k)*q[i]+(i==k)*q[j]))+ 
            5*(q[i]*Q(j,k)+q[j]*Q(i,k)+q[k]*Q(i,j));
          for(int l = k ; l < 3; ++l){
            // Octopole 
            H(i,j,k) += -2*q[l]*(Q(i,l)*(j==k)+Q(j,l)*(i==k)+Q(k,l)*(i==j));
            // Hexadecapole
            X(i,j,k,l) += m*(
              105.*q[i]*q[j]*q[k]*q[l]
              -15.*q2*(
                (i==j)*q[k]*q[l]+(i==l)*q[j]*q[k]+(i==k)*q[j]*q[l]+
                (j==l)*q[i]*q[k]+(j==k)*q[i]*q[l]+(l==k)*q[i]*q[j]
              )+3.*q4*((i==j)*(k==l)+(i==k)*(j==l)+(i==l)*(j==k)));
          }
        }
      }
    }
    return m; 
  }

  /**
  * @brief Compute the momenta for a center of mass in the tree. 
  * We gather two entities (particles or node) to start the process. 
  * Then we add remaining nodes and particles one by one moving the 
  * total mass and center of mass
  */
  void compute_momentum(
    tensor_u<double, symmetry_type::symmetric, 3, 3, 3, 3>& X,
    tensor_u<double, symmetry_type::symmetric, 3, 3, 3>& H,
    tensor_u<double, symmetry_type::symmetric, 3, 3>& Q,
    std::vector<body*> bs,
    std::vector<node*> ns)
  {
    int start_node = 0;
    int start_ent  = 0;
    point_t q{0}, c{0}; 
    double m, m0, m1, q2; 

    // Case of one sub-entity
    if(bs.size() + ns.size() == 1){
      if(bs.size() == 1){
        start_ent = 1; 
        q = bs[0]->coordinates();
        m0 = bs[0]->mass();
        m1 = 0;
      }else{
        start_node = 1; 
        q = ns[0]->coordinates();
        m0 = ns[0]->mass();
        m1 = 0;
        Q = ns[0]->quad(); 
        H = ns[0]->octo(); 
        X = ns[0]->hexa(); 
      }
    }else{
      if(bs.size() >= 2){
        start_ent = 2;
        // Combine two entities
        // This is also the parallel axis theorem
        q = bs[0]->coordinates()-bs[1]->coordinates();
        m0 = bs[0]->mass();
        m1 = bs[1]->mass();
        c = (m0*bs[0]->coordinates()+m1*bs[1]->coordinates())/(m0+m1);
      }else if(ns.size() >= 2){
        start_node = 2; 
        // Combine two nodes
        q = ns[0]->coordinates()-ns[1]->coordinates(); 
        m0 = ns[0]->mass(); 
        m1 = ns[1]->mass();
        Q = ns[0]->quad() + ns[1]->quad();  
        H = ns[0]->octo() + ns[1]->octo(); 
        X = ns[0]->hexa() + ns[1]->hexa(); 
        c = (m0*ns[0]->coordinates()+m1*ns[1]->coordinates())/(m0+m1);
      }else{
        start_node = 1; 
        start_ent = 1; 
        // Combine node and entity
        q = bs[0]->coordinates()-ns[0]->coordinates(); 
        m0 = bs[0]->mass(); 
        m1 = ns[0]->mass(); 
        Q = ns[0]->quad();
        H = ns[0]->octo();
        X = ns[0]->hexa(); 
        c = (m0*bs[0]->coordinates()+m1*ns[1]->coordinates())/(m0+m1); 
      }
    }
  
   //TODO : Remove this after everyone confirms that
   /*
    HL : I made a note here for adding expansions.

    1. Based on idea from ChanGA code
    First define shifting based on parallel axis theorem

    P_shift = m_tr (3*q^i*q^j - (i==j)*q^2)
    where m_tr = m1m2/(m1+m2) and same for H and X such that
    
    Then "adding" expansions will follow

    Q_tot = Q_1 + Q_2 + P_shift
    H_tot = H_1 + H_2 + P_shift
    X_tot = X_1 + X_2 + P_shift

    Note that there is index mismatching problem in this way

    2. Invent new way to add expansions
    It is advantageous to use composition formulae that incoporate
    with monopole, dipole, quadrupole, octopole, and hexadecapole.
    based on Benz et al. (1990). 
    P_shift from above is derived from composition formula for 
    quadrupole. For higher order, we may construct like

    (r_i - r_j) \tp (r_i - r_j) + (r_i - r_k) \tp (r_i - r_k) + (r_j - r_k) \tp (r_j - r_k)

    for octopole and where \tp is tensor product among different three position vector r^i

   */

    m = compute_XHQ(X,H,Q,m0,m1,q);
 
    // Loop over the entities remaining 
    for(int n = start_ent; n < bs.size(); ++n){
      m0 = m; 
      m1 = bs[n]->mass(); 
      q = bs[n]->coordinates()-c; 
      m = compute_XHQ(X,H,Q,m0,m1,q);
      c = (m0*c+m1*bs[n]->coordinates())/(m0+m1);
    }
    // Loop over the nodes remaining
    for(int n = start_node; n < ns.size(); ++n){
      m0 = m; 
      m1 = ns[n]->mass(); 
      q = ns[n]->coordinates()-c; 
      Q += ns[n]->quad(); 
      H += ns[n]->octo(); 
      X += ns[n]->hexa(); 
      m = compute_XHQ(X,H,Q,m0,m1,q);
      c = (m0*c+m1*ns[n]->coordinates())/(m0+m1);
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
    const point_t& local_coordinates,
    const point_t& dist_coordinates,
    const double& sm)
  {
    double dist = flecsi::distance(local_coordinates,dist_coordinates);
    point_t res = -gravitational_constant*sm/(dist*dist*dist)*
      (local_coordinates-dist_coordinates);
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
    const double& M, 
    const tensor_u<double, symmetry_type::symmetric, 3, 3, 3, 3>& X, 
    const tensor_u<double, symmetry_type::symmetric, 3, 3, 3>& H, 
    const tensor_u<double, symmetry_type::symmetric, 3, 3>& Q)
  {
    double d = flecsi::distance(local_coordinates,dist_coordinates);
    double d2 = d*d; 
    double d3 = d2*d; 
    double d5 = d3*d2; 
    double d7 = d5*d2; 
    double d9 = d7* d2; 
    double d11 = d9*d2; 
    point_t r = local_coordinates-dist_coordinates; 

    for(int m = 0; m < 3; ++m){
      // Monopole
      fc[m] += -gravitational_constant*M*r[m]/d3;
      #ifdef QUAD
      for(int i = 0 ; i < 3; ++i){
        // Quadrupole 
        fc[m] += Q(i,m)*r[i]/d5;
        for(int j = 0 ; j < 3; ++j){
          // Quadrupole 
          fc[m] += -2.5*Q(i,j)*r[i]*r[j]*r[m]/d7;
          #ifdef OCTO
          // Octopole 
          fc[m] += .5*H(i,j,m)*r[i]*r[j]/d7; 
          for(int k = 0 ; k < 3; ++k){
            // Octopole 
            fc[m] += -7./6.*H(i,j,k)*r[i]*r[j]*r[k]*r[m]/d9;
            #ifdef HEXA
            // Hexadecapole 
            fc[m] += 1./6.*X(i,j,k,m)*r[i]*r[j]*r[k]/d9;  
            for(int l = 0; l < 3; ++l){
              // Hexadecapole 
              fc[m] += -9./24.*X(i,j,k,l)*r[i]*r[j]*r[k]*r[l]*r[m]/d11;
            }
            #endif // HEXA
          }
          #endif // OCTO
        }
      }
      #endif // QUAD
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
    const double& M, 
    const tensor_u<double, symmetry_type::symmetric, 3, 3, 3, 3>& X, 
    const tensor_u<double, symmetry_type::symmetric, 3, 3, 3>& H, 
    const tensor_u<double, symmetry_type::symmetric, 3, 3>& Q)
  {
    double d = flecsi::distance(local_coordinates,dist_coordinates);
    double d2 = d*d;
    double d3 = d2*d;
    double d5 = d3*d2;
    double d7 = d5*d2;
    double d9 = d7* d2;
    double d11 = d9*d2;
    point_t r = local_coordinates-dist_coordinates; 

    for(int m = 0; m < 3; ++m){
      for(int q = 0; q < 3; ++q){
        int pr = m*3+q;
        // Monopole
        res[pr] += gravitational_constant*M/d3*(3.*r[m]*r[q]/d2-(q==m));
        #ifdef QUAD
        // Quadrupole
        res[pr] += Q(m,q)/d5; 
        for(int i = 0 ; i < 3; ++i){
          // Quadrupole
          res[pr] += -5.*(Q(i,m)*r[q]+Q(i,q)*r[m])*r[i]/d7; 
          #ifdef OCTO
          // Octopole 
          res[pr] += H(i,q,m)*r[i]/d7;
          #endif 
          for(int j = 0 ; j < 3; ++j){
            // Quadrupole 
            res[pr] += (35./2.*r[m]*r[q]/d2-2.5*(m==q))*
              Q(i,j)*r[i]*r[j]/d7;
            #ifdef OCTO 
            // Octopole
            res[pr] += -3.5*(H(i,j,m)*r[q]+H(i,j,q)*r[m])*r[i]*r[j]/d9;
            #ifdef HEXA
            // Hexadecapole
            res[pr] += .5*X(i,j,q,m)*r[i]*r[j]/d9; 
            #endif 
            for(int k = 0 ; k < 3; ++k){
              // Octopole 
              res[pr] += 7./6.*(9.*r[m]*r[q]/d2-(q==m))*
                H(i,j,k)*r[i]*r[j]*r[k]/d9; 
              #ifdef HEXA
              // Hexadecapole 
              res[pr] += -9./6.*(X(i,j,k,m)*r[q]+
                X(i,j,k,q)*r[m])*r[i]*r[j]*r[k]/d11;
              for(int l = 0; l < 3; ++l){
                // Hexadecapole 
                res[pr] += 9./24.*(11.*r[m]*r[q]/d2-(q==m))*
                  X(i,j,k,l)*r[i]*r[j]*r[k]*r[l]/d11;
              }
              #endif // HEXA
            }
            #endif // OCTO
          }
        }
        #endif // QUAD
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
    const double& M, 
    const tensor_u<double, symmetry_type::symmetric, 3, 3, 3, 3>& X, 
    const tensor_u<double, symmetry_type::symmetric, 3, 3, 3>& H, 
    const tensor_u<double, symmetry_type::symmetric, 3, 3>& Q)
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
    point_t r = local_coordinates-dist_coordinates; 

    for(int m = 0; m < 3; ++m){
      for(int q = 0; q < 3; ++q){
        for(int s = 0 ; s < 3; ++s){
          int pr = (m*3+q)*3+s;
          // Monopole 
          res[pr] += gravitational_constant*3.*M/d5*(
            (m==q)*r[s]+(q==s)*r[m]+(m==s)*r[q]-5.*r[m]*r[q]*r[s]/d2);
          #ifdef QUAD
          // Quadrupole
          res[pr] += -5./d7*
            (Q(m,q)*r[s]+Q(s,m)*r[q]+Q(s,q)*r[m]);
          #ifdef OCTO
          // Octopole 
          res[pr] += H(s,q,m)/d7;
          #endif
          for(int i = 0 ; i < 3; ++i){
            // Quadrupole 
            res[pr] += -(Q(i,m)*(q==s)+
                         Q(i,q)*(m==s)+
                         Q(i,s)*(m==q))*5.*r[i]/d7;
            res[pr] += (Q(i,m)*r[q]*r[s]+Q(i,q)*r[m]*r[s]+
              Q(i,s)*r[m]*r[q])*35.*r[i]/d9; 
            #ifdef OCTO
            // Octopole 
            res[pr] += -(H(i,q,m)*r[s]+
                         H(i,s,m)*r[q]+
                         H(i,s,q)*r[m])*7.*r[i]/d9;
            #endif 
            #ifdef HEXA
            // Hexadecapole 
            res[pr] += X(i,s,q,m)*r[i]/d9;
            #endif 
            for(int j = 0 ; j < 3; ++j){
              // Quadrupole
              res[pr] += 35./2.*((m==q)*r[s]+
                                 (m==s)*r[q]+
                                 (s==q)*r[m])*Q(i,j)*r[i]*r[j]/d9;
              res[pr] += -315./2.*Q(i,j)*r[i]*r[j]*r[m]*r[q]*r[s]/d11;
              #ifdef OCTO
              // Octopole 
              res[pr] += -3.5*(H(i,j,m)*(q==s)+
                               H(i,j,q)*(m==s)+
                               H(i,j,s)*(q==m))*r[i]*r[j]/d9;
              res[pr] += 63./2.*(H(i,j,m)*r[q]*r[s]+
                                 H(i,j,q)*r[m]*r[s]+
                                 H(i,j,s)*r[m]*r[q])*r[i]*r[j]/d11;
              #ifdef HEXA
              // Hexadecapole 
              res[pr] += -4.5*(X(i,j,q,m)*r[s]+
                               X(i,j,s,m)*r[q]+
                               X(i,j,s,q)*r[m])*r[i]*r[j]/d11;
              #endif 
              for(int k = 0 ; k < 3; ++k){
                // Octopole 
                res[pr] += 63./6.*((q==m)*r[s]+
                                   (m==s)*r[q]+
                                   (q==s)*r[m])*
                                   H(i,j,k)*r[i]*r[j]*r[k]/d11;
                res[pr] += -693./6.*H(i,j,k)*
                  r[i]*r[j]*r[k]*r[m]*r[q]*r[s]/d13;
                #ifdef HEXA
                // Hexadecapole 
                res[pr] += -9./6.*(X(i,j,k,m)*(q==s)+
                                   X(i,j,k,q)*(s==m)+
                                   X(i,j,k,s)*(q==m))*r[i]*r[j]*r[k]/d11;
                res[pr] += 99./6.*(X(i,j,k,m)*r[q]*r[s]+
                                   X(i,j,k,q)*r[m]*r[s]+
                                   X(i,j,k,s)*r[m]*r[q])*r[i]*r[j]*r[k]/d13;
                for(int l = 0; l < 3; ++l){
                  // Hexadecapole
                  res[pr] += 99./24.*((q==m)*r[s]+
                                      (m==s)*r[q]+
                                      (q==s)*r[m])*
                                      X(i,j,k,l)*r[i]*r[j]*r[k]*r[l]/d13;
                  res[pr] += -1287./24.*X(i,j,k,l)*
                      r[i]*r[j]*r[k]*r[l]*r[m]*r[q]*r[s]/d15; 
                }
                #endif // HEXA
              }
              #endif // OCTO
            }
          }
          #endif // QUAD
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


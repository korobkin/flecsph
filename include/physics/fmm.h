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

//#define OCTO
//#define HEXA

#pragma once

#include "params.h"
#include "tree.h"

namespace fmm {
  using namespace param;
  double gc = gravitational_constant;
  using sym_tensor_rank2 = flecsi::sym_tensor_rank2;
  using sym_tensor_rank3 = flecsi::sym_tensor_rank3;
  using sym_tensor_rank4 = flecsi::sym_tensor_rank4;

#if fmm_order == 1
  struct fmm_comms{
    using sym_tensor_rank2 = flecsi::sym_tensor_rank2;
    using sym_tensor_rank3 = flecsi::sym_tensor_rank3;
    using sym_tensor_rank4 = flecsi::sym_tensor_rank4;
    fmm_comms(){}
    fmm_comms(
      const point_t& _coords,
      const double& _T)
    : coords(_coords), T(_T){}
    fmm_comms( const node* n){
      coords = n->coordinates();
      T = n->mass();
    }
    void fill(const node* n){
      coords = n->coordinates();
      T = n->mass();
    }
    point_t coords;
    double T;
  };
#endif

#if fmm_order == 2
  struct fmm_comms {
    using sym_tensor_rank2 = flecsi::sym_tensor_rank2;
    using sym_tensor_rank3 = flecsi::sym_tensor_rank3;
    using sym_tensor_rank4 = flecsi::sym_tensor_rank4;
    fmm_comms(){}
    fmm_comms(
      const point_t& _coords,
      const double& _T,
      const sym_tensor_rank2& _Q)
    : coords(_coords), T(_T), Q(_Q){}
    fmm_comms( const node* n){
      coords = n->coordinates();
      T = n->mass();
      Q = n->quad();
    }
    void fill(const node* n){
      coords = n->coordinates();
      T = n->mass();
      Q = n->quad();
    }
    point_t coords;
    double T;
    sym_tensor_rank2 Q;
  };
#endif // fmm_order == 2

#if fmm_order == 3
  struct fmm_comms {
    using sym_tensor_rank2 = flecsi::sym_tensor_rank2;
    using sym_tensor_rank3 = flecsi::sym_tensor_rank3;
    using sym_tensor_rank4 = flecsi::sym_tensor_rank4;
    fmm_comms(){}
    fmm_comms(
      const point_t& _coords,
      const double& _T,
      const sym_tensor_rank3& _H,
      const sym_tensor_rank2& _Q)
    : coords(_coords), T(_T), H(_H), Q(_Q){}
    fmm_comms( const node* n){
      coords = n->coordinates();
      T = n->mass();
      H = n->octo();
      Q = n->quad();
    }
    void fill(const node* n){
      coords = n->coordinates();
      T = n->mass();
      H = n->octo();
      Q = n->quad();
    }
    point_t coords;
    double T;
    sym_tensor_rank3 H;
    sym_tensor_rank2 Q;
  };
#endif // fmm_order == 3

#if fmm_order == 4
  struct fmm_comms {
    using sym_tensor_rank2 = flecsi::sym_tensor_rank2;
    using sym_tensor_rank3 = flecsi::sym_tensor_rank3;
    using sym_tensor_rank4 = flecsi::sym_tensor_rank4;
    fmm_comms(){}
    fmm_comms(
      const point_t& _coords,
      const double& _T,
      const sym_tensor_rank4& _X,
      const sym_tensor_rank3& _H,
      const sym_tensor_rank2& _Q)
    : coords(_coords), T(_T), X(_X), H(_H), Q(_Q){}
    fmm_comms( const node* n){
      coords = n->coordinates();
      T = n->mass();
      X = n->hexa();
      H = n->octo();
      Q = n->quad();
    }
    void fill(const node* n){
      coords = n->coordinates();
      T = n->mass();
      X = n->hexa();
      H = n->octo();
      Q = n->quad();
    }
    point_t coords;
    double T;
    sym_tensor_rank4 X;
    sym_tensor_rank3 H;
    sym_tensor_rank2 Q;
  };
#endif



  /**
  * Sum quadrupole moment for two bodies
  */
  inline double compute_Q(
    // Left moments = where we sum
    sym_tensor_rank2& Ql,
    // Right moments
    const sym_tensor_rank2& Qr,
    // Left and right masses
    const double& ml, const double& mr,
    // Left and right positions
    const point_t& pl,const point_t& pr)
  {
    // q = left - right
    const point_t q = pl-pr;
    const double q2 = q[0]*q[0]+q[1]*q[1]+q[2]*q[2];
    const double q4 = q2*q2;
    // Reduced mass and moments
    const sym_tensor_rank2 r_Q = (mr*Ql-ml*Qr)/(ml+mr);
    const double r_m = ml*mr/(ml+mr);
    // We sum the result on left: Ql, Hl and Xl
    for(int i = 0 ; i < 3; ++i){
      for(int j = i ; j < 3; ++j){
        // Quadrupole
        Ql(i,j) += r_m*(3.*q[i]*q[j] - (i==j)*q2);
      }
    }
    // Add the right component to left
    Ql += Qr;
    return r_m;
  }

  /**
  * Sum octupole and quadrupole moments for two bodies
  */
  inline double compute_HQ(
    // Left moments = where we sum
    sym_tensor_rank3& Hl,
    sym_tensor_rank2& Ql,
    // Right moments
    const sym_tensor_rank3& Hr,
    const sym_tensor_rank2& Qr,
    // Left and right masses
    const double& ml, const double& mr,
    // Left and right positions
    const point_t& pl,const point_t& pr)
  {
    // q = left - right
    const point_t q = pl - pr;
    const double q2 = q[0]*q[0]+q[1]*q[1]+q[2]*q[2];
    const double q4 = q2*q2;
    // Reduced mass and moments
    const sym_tensor_rank2 r_Q = (mr*Ql - ml*Qr)/(ml + mr);
    const double r_m = ml*mr/(ml + mr);
    // We sum the result on left: Ql, Hl and Xl
    for(int i = 0 ; i < 3; ++i){
      for(int j = i ; j < 3; ++j){
        // Quadrupole
        Ql(i,j) += r_m*(3.*q[i]*q[j]-(i==j)*q2);
        for(int k = j ; k < 3; ++k){
          // Octopole
          Hl(i,j,k) += r_m*((mr-ml)/(mr+ml))*(
            15.*q[i]*q[j]*q[k]
            -3.*q2*((i==j)*q[k]+(j==k)*q[i]+(i==k)*q[j]))+
            5.*(q[i]*r_Q(j,k)+q[j]*r_Q(i,k)+q[k]*r_Q(i,j));
          for(int s = 0 ; s < 3; ++s){
            Hl(i,j,k) += -2*q[s]*(r_Q(i,s)*(j==k)+r_Q(j,s)*(i==k)+r_Q(k,s)*(i==j));
          }
        }
      }
    }
    // Add the right component to left
    Hl += Hr;
    Ql += Qr;

    #if 0
    for(int l = k ; l < 3; ++l){
      // Hexadecapole
      X(i,j,k,l) += m*(
        105.*q[i]*q[j]*q[k]*q[l]
        -15.*q2*(
          (i==j)*q[k]*q[l]+(i==l)*q[j]*q[k]+(i==k)*q[j]*q[l]+
          (j==l)*q[i]*q[k]+(j==k)*q[i]*q[l]+(l==k)*q[i]*q[j]
        )+3.*q4*((i==j)*(k==l)+(i==k)*(j==l)+(i==l)*(j==k)));
    }
    #endif
    return r_m;
  }


  /**
  * Sum quad, octo and hexapoles for two entities
  */
  inline double compute_XHQ(
    // Left moments = where we sum
    sym_tensor_rank4& Xl,
    sym_tensor_rank3& Hl,
    sym_tensor_rank2& Ql,
    // Right moments
    const sym_tensor_rank4& Xr,
    const sym_tensor_rank3& Hr,
    const sym_tensor_rank2& Qr,
    // Left and right masses
    const double& ml, const double& mr,
    // Left and right positions
    const point_t& pl,const point_t& pr)
  {
    // q = left - right
    const point_t q = pl-pr;
    const double q2 = q[0]*q[0]+q[1]*q[1]+q[2]*q[2];
    const double q4 = q2*q2;
    // Reduced mass and moments
    const sym_tensor_rank2 r_Q = (mr*Ql-ml*Qr)/(ml+mr);
    const double r_m = ml*mr/(ml+mr);
    // We sum the result on left: Ql, Hl and Xl
    for(int i = 0 ; i < 3; ++i){
      for(int j = i ; j < 3; ++j){
        // Quadrupole
        Ql(i,j) += r_m*(3.*q[i]*q[j]-(i==j)*q2);
        for(int k = j ; k < 3; ++k){
          // Octopole
          Hl(i,j,k) += r_m*((mr-ml)/(mr+ml))*(
            15.*q[i]*q[j]*q[k]
            -3.*q2*((i==j)*q[k]+(j==k)*q[i]+(i==k)*q[j]))+
            5.*(q[i]*r_Q(j,k)+q[j]*r_Q(i,k)+q[k]*r_Q(i,j));
          for(int s = 0 ; s < 3; ++s){
            Hl(i,j,k) += -2*q[s]*(r_Q(i,s)*(j==k)+r_Q(j,s)*(i==k)+r_Q(k,s)*(i==j));
          }
        }
      }
    }
    // Add the right component to left
    Xl += Xr;
    Hl += Hr;
    Ql += Qr;
    // Xl, Ql and Hl now contains total moment left + right

    #if 0
    for(int l = k ; l < 3; ++l){
      // Hexadecapole
      X(i,j,k,l) += m*(
        105.*q[i]*q[j]*q[k]*q[l]
        -15.*q2*(
          (i==j)*q[k]*q[l]+(i==l)*q[j]*q[k]+(i==k)*q[j]*q[l]+
          (j==l)*q[i]*q[k]+(j==k)*q[i]*q[l]+(l==k)*q[i]*q[j]
        )+3.*q4*((i==j)*(k==l)+(i==k)*(j==l)+(i==l)*(j==k)));
    }
    #endif
    return r_m;
  }


  /**
  * @brief Compute the moments for a center of mass in the tree.
  *        We gather two entities (particles or nodes) to start the process.
  *        Then we add remaining nodes and particles one by one moving the
  *        total mass and center of mass
  */
  void compute_moments(
    node* cofm,
    std::vector<body*>& bs,
    std::vector<node*>& ns)
  {

    #if fmm_order > 1 // if fmm_order is 1, this function is empty
    sym_tensor_rank2& Q = cofm->quad();
    #if fmm_order > 2
      sym_tensor_rank3& H = cofm->octo();
      #if fmm_order > 3
        sym_tensor_rank4& X = cofm->hexa();
      #endif
    #endif

    // zero out Taylor expansion coefficients
    cofm->pc() = 0;
    cofm->fc() = 0;
    cofm->dfcdr() = 0;
    #if fmm_order > 2
      cofm->dfcdrdr() = 0;
      #if fmm_order > 3
        cofm->dfcdrdrdr() = 0;
      #endif
    #endif

    int start_node = 0;
    int start_ent  = 0;
    // COM informations
    point_t c{0};
    double m;

    // Start from a particle
    // Make this particle the current COM
    if(bs.size() > 0){
      start_ent = 1;
      c = bs[0]->coordinates();
      m = bs[0]->mass();
      Q = {0};
      #if fmm_order > 2
        H = {0};
        #if fmm_order > 3
          X = {0};
        #endif
      #endif
    }else{
      start_node = 1;
      // Start from a node
      // Make this node the current COM
      c = ns[0]->coordinates();
      m = ns[0]->mass();
      Q = ns[0]->quad();
      #if fmm_order > 2
        H = ns[0]->octo();
        #if fmm_order > 3
          X = ns[0]->hexa();
        #endif
      #endif
    }

    // For particles: moments are = 0
    const sym_tensor_rank4 Xp = {0};
    const sym_tensor_rank3 Hp = {0};
    const sym_tensor_rank2 Qp = {0};
    // Add particles
    for(int k = start_ent; k < bs.size(); ++k){
      const body & b = *(bs[k]);
      const double mp = b.mass();
      const point_t & cp = b.coordinates();
      #if fmm_order == 2
        compute_Q(Q, Qp, m, mp, c, cp);
      #elif fmm_order == 3
        compute_HQ(H, Q, Hp, Qp, m, mp, c, cp);
      #elif fmm_order == 4
        compute_XHQ(X, H, Q, Xp, Hp, Qp, m, mp, c, cp);
      #else
        assert(false);
      #endif
      // New COM mass and coordinates
      c = (m*c + mp*cp)/(m + mp);
      m += mp;
    }
    // Add nodes
    for(int k = start_node; k < ns.size(); ++k){
      const node & n = *(ns[k]);
      const double mn = n.mass();
      const point_t & cn = n.coordinates();
      #if fmm_order == 2
        compute_Q(Q, n.quad(), m, mn, c, cn);
      #elif fmm_order == 3
        compute_HQ(H, Q, n.octo(), n.quad(), m, mn, c, cn);
      #elif fmm_order == 4
        compute_XHQ(X, H, Q, n.hexa(), n.octo(), n.quad(), m, mn, c, cn);
      #else
        assert(false);
      #endif
      // New COM mass and coordinates
      c = (m*c + mn*cn)/(m + mn);
      m += mn;
    }
    #endif // if constexpr gdimension == 3, fmm_order > 1
  }


  /*
  * @brief Compute the gravitation interaction
  * Return the computed value if needed for direct particle interaction
  *
  * The Sink is the one on which I compute the fc, dfcdr, dfcdrdr, dfcdrdrdr
  * The Source is the distant box
  */
  inline
  point_t gravitation_p2p(
    double& gpot,
    const point_t& local_coordinates,
    const point_t& dist_coordinates,
    const double& sm)
  {
    double dist = flecsi::distance(local_coordinates,dist_coordinates);
    gpot += -gc*sm/dist;
    point_t res = -gc*sm/(dist*dist*dist)*
      (local_coordinates-dist_coordinates);
    return res;
  }


  /*
  * @brief Compute the gravitation interaction
  * Return the computed value if needed for direct particle interaction
  *
  * The Sink is the one on which I compute the fc, dfcdr, dfcdrdr
  * The Source is the distant box
  */
  inline
  void gravitation_fc(
    double& pc,
    point_t& fc,
    const point_t& local_coordinates,
    const fmm_comms& source)
  {
    const point_t& dist_coordinates = source.coords;
    const double M = source.T;
    #if fmm_order > 1
      const sym_tensor_rank2& Q = source.Q;
      #if fmm_order > 2
        const sym_tensor_rank3& H = source.H;
        #if fmm_order > 3
          const sym_tensor_rank4& X = source.X;
        #endif
      #endif
    #endif

    double d = flecsi::distance(local_coordinates,dist_coordinates);
    double d2 = d*d;
    double d3 = d2*d;
    double d5 = d3*d2;
    double d7 = d5*d2;
    double d9 = d7* d2;
    double d11 = d9*d2;
    point_t r = local_coordinates-dist_coordinates;

    pc += -gc*M/d;

    for(int m = 0; m < gdimension; ++m){
      // Monopole
      fc[m] += -gc*M*r[m]/d3;
      #if fmm_order > 1
      for(int i = 0 ; i < gdimension; ++i){
        // Quadrupole contribution to the potential
        pc += -gc*.5*Q(m,i)*r[m]*r[i]/d5;

        #if fmm_order > 2
        // Quadrupole
        fc[m] += gc*Q(i,m)*r[i]/d5;
        for(int j = 0 ; j < gdimension; ++j){
          // Quadrupole
          fc[m] += -gc*2.5*Q(i,j)*r[i]*r[j]*r[m]/d7;

          // Octopole Potential
          pc += -gc*1./6.*H(m,i,j)*r[m]*r[i]*r[j]/d7;
          #if fmm_order > 3
          // Octopole
          fc[m] += gc*.5*H(i,j,m)*r[i]*r[j]/d7;
          for(int k = 0 ; k < gdimension; ++k){
            // Octopole
            fc[m] += -gc*7./6.*H(i,j,k)*r[i]*r[j]*r[k]*r[m]/d9;
            // Hexadecapole Potential
            pc += -gc*1./24.*X(m,i,j,k)*r[m]*r[i]*r[j]*r[k]/d9;
            #if fmm_order > 4
            // Hexadecapole
            fc[m] += gc*1./6.*X(i,j,k,m)*r[i]*r[j]*r[k]/d9;
            for(int l = 0; l < gdimension; ++l){
              // Hexadecapole
              fc[m] += -gc*9./24.*X(i,j,k,l)*r[i]*r[j]*r[k]*r[l]*r[m]/d11;
            }
            #endif // fmm_order > 4
          }
          #endif // fmm_order > 3
        }
        #endif // fmm_order > 2
      } // for i
      #endif // fmm_order > 1
    } // for m
  }


  /*
  * @brief Compute the Jacobian (dfcdr) matrix on the Sink from the Source
  */
  inline
  void gravitation_dfcdr(
    sym_tensor_rank2& res,
    const point_t& local_coordinates,
    const fmm_comms& source)
  {

    if constexpr (fmm_order <= 1) return;
    const point_t& dist_coordinates = source.coords;
    const double M = source.T;
    #if fmm_order > 1
      const sym_tensor_rank2& Q = source.Q;
      #if fmm_order > 2
        const sym_tensor_rank3& H = source.H;
        #if fmm_order > 3
          const sym_tensor_rank4& X = source.X;
        #endif
      #endif
    #endif

    double d = flecsi::distance(local_coordinates,dist_coordinates);
    double d2 = d*d;
    double d3 = d2*d;
    double d5 = d3*d2;
    double d7 = d5*d2;
    double d9 = d7* d2;
    double d11 = d9*d2;
    point_t r = local_coordinates-dist_coordinates;

    for(int m = 0; m < gdimension; ++m){
      for(int q = m; q < gdimension; ++q){
        // Monopole
        res(m,q) +=gc* M/d3*(3.*r[m]*r[q]/d2-(q==m));

        #if fmm_order > 2
        // Quadrupole
        res(m,q) += gc*Q(m,q)/d5;
        for(int i = 0 ; i < gdimension; ++i){
          // Quadrupole
          res(m,q) += -gc*5.*(Q(i,m)*r[q]+Q(i,q)*r[m])*r[i]/d7;
          for(int j = 0 ; j < gdimension; ++j){
            // Quadrupole
            res(m,q) += gc*(35./2.*r[m]*r[q]/d2-2.5*(m==q))*
              Q(i,j)*r[i]*r[j]/d7;
          }

          #if fmm_order > 3
          // Octopole
          res(m,q) += gc*H(i,q,m)*r[i]/d7;
          for(int j = 0 ; j < gdimension; ++j){
            // Octopole
            res(m,q) += -gc*3.5*(H(i,j,m)*r[q]+H(i,j,q)*r[m])*r[i]*r[j]/d9;
            #if fmm_order > 4
            // Hexadecapole
            res(m,q) += gc*.5*X(i,j,q,m)*r[i]*r[j]/d9;
            #endif
            for(int k = 0 ; k < gdimension; ++k){
              // Octopole
              res(m,q) += gc*7./6.*(9.*r[m]*r[q]/d2-(q==m))*
                H(i,j,k)*r[i]*r[j]*r[k]/d9;
              #if fmm_order > 4
              // Hexadecapole
              res(m,q) += -gc*9./6.*(X(i,j,k,m)*r[q]+
                X(i,j,k,q)*r[m])*r[i]*r[j]*r[k]/d11;
              for(int l = 0; l < gdimension; ++l){
                // Hexadecapole
                res(m,q) += gc*9./24.*(11.*r[m]*r[q]/d2-(q==m))*
                  X(i,j,k,l)*r[i]*r[j]*r[k]*r[l]/d11;
              }
              #endif // fmm_order > 4
            } // for k
          } // for j
          #endif // fmm_order > 3
        } // for i
        #endif // fmm_order > 2
      } // for q
    } // for m
  }

  /*
  * @brief Compute the Hessian (dfcdrdr) matrix on the Sink from the Source
  */
  inline
  void gravitation_dfcdrdr(
    sym_tensor_rank3& res,
    const point_t& local_coordinates,
    const fmm_comms& source)
  {

    if constexpr (fmm_order <= 2) return;
    const point_t& dist_coordinates = source.coords;
    const double M = source.T;
    #if fmm_order > 1
      const sym_tensor_rank2& Q = source.Q;
      #if fmm_order > 2
        const sym_tensor_rank3& H = source.H;
        #if fmm_order > 3
          const sym_tensor_rank4& X = source.X;
        #endif
      #endif
    #endif

    double d = flecsi::distance(local_coordinates,dist_coordinates);
    point_t r = local_coordinates-dist_coordinates;
    const double d2 = d*d;
    const double d3 = d2*d;
    const double d5 = d3*d2;
    #if fmm_order == 3
      for(int m = 0; m < gdimension; ++m)
      for(int q = m; q < gdimension; ++q)
      for(int s = q; s < gdimension; ++s) {
        res(m,q,s) += gc*3.*M/d5*(
          (m==q)*r[s]+(q==s)*r[m]+(m==s)*r[q]-5.*r[m]*r[q]*r[s]/d2);
      }
    #endif

    #if fmm_order > 3
    const double d7 = d5*d2;
    const double d9 = d7*d2;
    const double d11 = d9*d2;
    const double d13 = d11*d2;
    const double d15 = d13*d2;

    for(int m = 0; m < gdimension; ++m){
      for(int q = m; q < gdimension; ++q){
        for(int s = q ; s < gdimension; ++s){
          // Monopole
          res(m,q,s) += gc*3.*M/d5*(
            (m==q)*r[s]+(q==s)*r[m]+(m==s)*r[q]-5.*r[m]*r[q]*r[s]/d2);
          // Quadrupole
          res(m,q,s) += -gc*5./d7*
            (Q(m,q)*r[s]+Q(s,m)*r[q]+Q(s,q)*r[m]);
          // Octopole
          res(m,q,s) += gc*H(s,q,m)/d7;
          for(int i = 0 ; i < gdimension; ++i){
            // Quadrupole
            res(m,q,s) += -gc*(Q(i,m)*(q==s)+
                         Q(i,q)*(m==s)+
                         Q(i,s)*(m==q))*5.*r[i]/d7;
            res(m,q,s) += gc*(Q(i,m)*r[q]*r[s]+Q(i,q)*r[m]*r[s]+
              Q(i,s)*r[m]*r[q])*35.*r[i]/d9;
            // Octopole
            res(m,q,s) += -gc*(H(i,q,m)*r[s]+
                         H(i,s,m)*r[q]+
                         H(i,s,q)*r[m])*7.*r[i]/d9;
            #ifdef HEXA
            // Hexadecapole
            res(m,q,s) += gc*X(i,s,q,m)*r[i]/d9;
            #endif
            for(int j = 0 ; j < gdimension; ++j){
              // Quadrupole
              res(m,q,s) += gc*35./2.*((m==q)*r[s]+
                                 (m==s)*r[q]+
                                 (s==q)*r[m])*Q(i,j)*r[i]*r[j]/d9;
              res(m,q,s) += -gc*315./2.*Q(i,j)*r[i]*r[j]*r[m]*r[q]*r[s]/d11;
              // Octopole
              res(m,q,s) += -gc*3.5*(H(i,j,m)*(q==s)+
                               H(i,j,q)*(m==s)+
                               H(i,j,s)*(q==m))*r[i]*r[j]/d9;
              res(m,q,s) += gc*63./2.*(H(i,j,m)*r[q]*r[s]+
                                 H(i,j,q)*r[m]*r[s]+
                                 H(i,j,s)*r[m]*r[q])*r[i]*r[j]/d11;
              #ifdef HEXA
              // Hexadecapole
              res(m,q,s) += -gc*4.5*(X(i,j,q,m)*r[s]+
                               X(i,j,s,m)*r[q]+
                               X(i,j,s,q)*r[m])*r[i]*r[j]/d11;
              #endif
              for(int k = 0 ; k < gdimension; ++k){
                // Octopole
                res(m,q,s) += gc*63./6.*((q==m)*r[s]+
                                   (m==s)*r[q]+
                                   (q==s)*r[m])*
                                   H(i,j,k)*r[i]*r[j]*r[k]/d11;
                res(m,q,s) += -gc*693./6.*H(i,j,k)*
                  r[i]*r[j]*r[k]*r[m]*r[q]*r[s]/d13;
                #ifdef HEXA
                // Hexadecapole
                res(m,q,s) += -gc*9./6.*(X(i,j,k,m)*(q==s)+
                                   X(i,j,k,q)*(s==m)+
                                   X(i,j,k,s)*(q==m))*r[i]*r[j]*r[k]/d11;
                res(m,q,s) += gc*99./6.*(X(i,j,k,m)*r[q]*r[s]+
                                   X(i,j,k,q)*r[m]*r[s]+
                                   X(i,j,k,s)*r[m]*r[q])*r[i]*r[j]*r[k]/d13;
                for(int l = 0; l < gdimension; ++l){
                  // Hexadecapole
                  res(m,q,s) += gc*99./24.*((q==m)*r[s]+
                                      (m==s)*r[q]+
                                      (q==s)*r[m])*
                                      X(i,j,k,l)*r[i]*r[j]*r[k]*r[l]/d13;
                  res(m,q,s) += -gc*1287./24.*X(i,j,k,l)*
                      r[i]*r[j]*r[k]*r[l]*r[m]*r[q]*r[s]/d15;
                }
                #endif // HEXA
              }
            } // for j
          } // for i
        } // for s
      } // for q
    } // for m
    #endif // fmm_order > 3
  } //dfcdrdr


  /*
  * @brief Compute the D^3 f (dfcdrdrdr) matrix on the Sink from the Source
  *  TODO : We only need monopole for this but need to check
  */
  inline
  void gravitation_dfcdrdrdr(
    sym_tensor_rank4& res,
    const point_t& local_coordinates,
    const fmm_comms& source)
  {
    #if fmm_order > 3
    const point_t& dist_coordinates = source.coords;
    const double M = source.T;
    double d = flecsi::distance(local_coordinates,dist_coordinates);
    const double d2 = d*d;
    const double d4 = d2*d2;
    const double d5 = d4*d;
    point_t r = local_coordinates-dist_coordinates;

    for(int i = 0; i < gdimension; ++i){
      for(int j = i; j < gdimension; ++j){
        for(int k = j; k < gdimension; ++k){
         for(int l = k; l < gdimension; ++l){
          // Monopole
          res(i,j,k,l) += gc*3.*M/d5*(
            (i==j)*(k==l)+(j==k)*(i==l)+(i==k)*(j==l)
            -5./d2*((i==j)*r[k]*r[l]+(j==k)*r[i]*r[l]+(i==k)*r[j]*r[l]+
                    (i==l)*r[j]*r[k]+(j==l)*r[i]*r[k]+(k==l)*r[i]*r[j])
            +35./d4*r[i]*r[j]*r[k]*r[l]);
              } // for l
            } // for k
          } // for j
    } // for i
    #endif // fmm_order > 3

  } // dfcdrdrdr

  /*
  * @brief Taylor expansion of degree 2 using the computed function, jacobi,
  * hessian and the targeted particle
  */
  void interaction_c2p(
    body* sink,
    const node* source)
  {
    const double& pc  = source->pc();
    const point_t& fc = source->fc();
    point_t cofm_coordinates = source->coordinates();

    point_t part_coordinates = sink->coordinates();
    point_t r = part_coordinates - cofm_coordinates;
    point_t grav = fc;
    double pot = pc;

    for(int i = 0 ; i < gdimension; ++i){
      pot += -r[i]*fc[i];
    }

    #if fmm_order > 1
      const sym_tensor_rank2& dfcdr = source->dfcdr();
      // The Jacobi
      for(int m=0;m<gdimension;++m){
        for(int i=0;i<gdimension;++i){
          grav[m] += dfcdr(m,i)*r[i];
          pot += -.5*r[m]*r[i]*dfcdr(m,i);
        } // for i
      } // for m
    #endif

    #if fmm_order > 2
      const sym_tensor_rank3& dfcdrdr = source->dfcdrdr();
      // The hessian
      for(int m = 0; m < gdimension; ++m){
        for(int i = 0; i < gdimension; ++i){
          for(int j = 0 ; j < gdimension; ++j){
            grav[m] += .5*r[i]*r[j]*dfcdrdr(m,i,j);
            pot += -1./6.*r[m]*r[i]*r[j]*dfcdrdr(m,i,j);
          } // for j
        } // for i
      } // for m
    #endif

    #if fmm_order > 3
      const sym_tensor_rank4& dfcdrdrdr = source->dfcdrdrdr();
      // The D^3 f
      for(int m = 0; m < gdimension; ++m){
        for(int i = 0; i < gdimension; ++i){
          for(int j = 0; j < gdimension; ++j){
           for(int k = 0; k < gdimension; ++k){
             grav[m] += 1./.5*r[i]*r[j]*r[k]*dfcdrdrdr(m,i,j,k);
             pot += -1./24.*r[m]*r[i]*r[j]*r[k]*dfcdrdrdr(m,i,j,k);
            } // for k
          } // for j
        } // for i
      } // for m
    #endif

    sink->setGPotential(sink->getGPotential()+pot);
    sink->setGAcceleration(grav+sink->getGAcceleration());
  }


  /**
  * @brief For all sub_entities, add gravitational force and potential of the
  *        node n, using Taylor expansion coefficients stored in the node.
  */
  void fmm_c2p(
    const node* nd,
    std::vector<body*>& sub_entities)
  {
    for (int k = 0; k < sub_entities.size(); ++k) {
      interaction_c2p(sub_entities[k], nd);
    } // for
  } // fmm_c2p


  /**
  * @brief Particle-particle interactions between 'sources' and 'sinks'
  */
  void fmm_p2p(
    std::vector<body*> &sinks,
    const std::vector<fmm_comms>& node_sources,
    const std::vector<body*>& particle_sources)
  {
    for (int i=0; i<sinks.size(); ++i) {
      body *p = sinks[i];
      double pc = p->getGPotential();
      point_t acc = p->getGAcceleration();
      for (int k = 0; k < node_sources.size(); ++k) {
        gravitation_fc(pc, acc, p->coordinates(), node_sources[k]);
      } // for
      for (int k = 0; k < particle_sources.size(); ++k) {
        body *q = particle_sources[k];
        if (q->id() == p->id())
          continue;
        acc+= gravitation_p2p(pc,p->coordinates(),q->coordinates(),q->mass());
      } // for
      p->setGPotential(pc);
      p->setGAcceleration(acc);
    }
  }


  /**
  * @brief node-node interaction: update Taylor expansion coefficients
  */
  void taylor_c2c(
    node* sink,
    const node* source)
  {
    fmm_comms fmm_source{source};
    gravitation_fc(sink->pc(), sink->fc(), sink->coordinates(), fmm_source);
    #if fmm_order > 1
      gravitation_dfcdr(sink->dfcdr(), sink->coordinates(), fmm_source);
      #if fmm_order > 2
        gravitation_dfcdrdr(sink->dfcdrdr(), sink->coordinates(), fmm_source);
      #endif
    #endif
  }

  /**
  * @brief node<-particle interaction: update Taylor expansion coefficients
  */
  void taylor_p2c(
    node* sink,
    const type_t point_source_mass,
    const point_t& point_source_coords)
  {
    fmm_comms source;
    source.coords = point_source_coords;
    source.T = point_source_mass;
    #if fmm_order > 1
      source.Q = 0;
      #if fmm_order > 2
        source.H = 0;
        #if fmm_order > 3
          source.X = 0;
        #endif
      #endif
    #endif

    gravitation_fc(sink->pc(), sink->fc(), sink->coordinates(), source);
    #if fmm_order > 1
      gravitation_dfcdr(sink->dfcdr(), sink->coordinates(), source);
      #if fmm_order > 2
        gravitation_dfcdrdr(sink->dfcdrdr(), sink->coordinates(), source);
        #if fmm_order > 3
          gravitation_dfcdrdrdr(sink->dfcdrdrdr(), sink->coordinates(), source);
        #endif
      #endif
    #endif
  }

} // namespace fmm


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

  struct fmm_comms{
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


  /**
  * Compute momenta for two entities
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
  * @brief Compute the momenta for a center of mass in the tree. 
  * We gather two entities (particles or node) to start the process. 
  * Then we add remaining nodes and particles one by one moving the 
  * total mass and center of mass
  */
  void compute_moments(
    node* cofm,
    std::vector<body*>& bs,
    std::vector<node*>& ns)
  {

    if constexpr (gdimension == 3 and fmm_order > 1) {
      sym_tensor_rank4& X = cofm->hexa();
      sym_tensor_rank3& H = cofm->octo();
      sym_tensor_rank2& Q = cofm->quad();

      // zero out Taylor expansion coefficients
      cofm->pc() = 0;
      cofm->fc() = 0;
      cofm->dfcdr() = 0;
      cofm->dfcdrdr() = 0;
      
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
        H = {0}; 
        X = {0}; 
      }else{
        start_node = 1; 
        // Start from a node 
        // Make this node the current COM 
        c = ns[0]->coordinates();
        m = ns[0]->mass(); 
        Q = ns[0]->quad();  
        H = ns[0]->octo(); 
        X = ns[0]->hexa();  
      }

      // For particles: moments are = 0 
      const sym_tensor_rank4 Xp = {0};
      const sym_tensor_rank3 Hp = {0}; 
      const sym_tensor_rank2 Qp = {0};
      // Add particles 
      for(int n = start_ent; n < bs.size(); ++n){
        compute_XHQ(
          // Left moments = where we sum 
          X,H,Q,
          // Rights moments = particles = 0 
          Xp,Hp,Qp,
          // Masses
          m,bs[n]->mass(),
          // Coordinates
          c,bs[n]->coordinates());
        // New COM mass and coordinates 
        c = (m*c+bs[n]->mass()*bs[n]->coordinates())/(m+bs[n]->mass());
        m += bs[n]->mass(); 
      }
      // Add nodes
      for(int n = start_node; n < ns.size(); ++n){
        compute_XHQ(
          // Left moments = where we sum 
          X,H,Q,
          // Right moments = this node moments
          ns[n]->hexa(),ns[n]->octo(),ns[n]->quad(),
          // Masses
          m,ns[n]->mass(),
          // Coordinates 
          c,ns[n]->coordinates());
        // New COM mass and coordinates 
        c = (m*c+ns[n]->mass()*ns[n]->coordinates())/(m+ns[n]->mass());
        m += ns[n]->mass();
      }
    } // if constexpr gdimension == 3, fmm_order > 1
  }


  /*
  * @brief Compute the gravitation interaction
  * Return the computed value if needed for direct particle interaction
  *
  * The Sink is the one on which I compute the fc, dfcdr, dfcdrdr
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
    const point_t& dist_coordinates,
    const double& M, 
    const sym_tensor_rank4& X, 
    const sym_tensor_rank3& H, 
    const sym_tensor_rank2& Q)
  {
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
      if constexpr(fmm_order <= 1) continue;
      
      for(int i = 0 ; i < gdimension; ++i){
        // Quadrupole contribution to the potential 
        pc += -gc*.5*Q(m,i)*r[m]*r[i]/d5;
        if constexpr(fmm_order == 2) continue;
        
        // Quadrupole 
        fc[m] += gc*Q(i,m)*r[i]/d5;
        for(int j = 0 ; j < gdimension; ++j){
          // Quadrupole 
          fc[m] += -gc*2.5*Q(i,j)*r[i]*r[j]*r[m]/d7;

          // Octopole Potential 
          pc += -gc*1./6.*H(m,i,j)*r[m]*r[i]*r[j]/d7;
          if constexpr(fmm_order == 3) continue;
          // Octopole 
          fc[m] += gc*.5*H(i,j,m)*r[i]*r[j]/d7; 
          for(int k = 0 ; k < gdimension; ++k){
            // Octopole 
            fc[m] += -gc*7./6.*H(i,j,k)*r[i]*r[j]*r[k]*r[m]/d9;
            #ifdef HEXA
            // Hexadecapole Potential 
            pc += -gc*1./24.*X(m,i,j,k)*r[m]*r[i]*r[j]*r[k]/d9; 
            // Hexadecapole 
            fc[m] += gc*1./6.*X(i,j,k,m)*r[i]*r[j]*r[k]/d9;  
            for(int l = 0; l < gdimension; ++l){
              // Hexadecapole 
              fc[m] += -gc*9./24.*X(i,j,k,l)*r[i]*r[j]*r[k]*r[l]*r[m]/d11;
            }
            #endif // HEXA
          }
        }
      } // for i
    } // for m
  }


  /*
  * @brief Compute the Jacobian (dfcdr) matrix on the Sink from the Source
  */
  inline
  void gravitation_dfcdr(
    sym_tensor_rank2& res,
    const point_t& local_coordinates,
    const point_t& dist_coordinates,
    const double& M, 
    const sym_tensor_rank4& X, 
    const sym_tensor_rank3& H, 
    const sym_tensor_rank2& Q)
  {
    
    if constexpr (fmm_order <= 1) return;
    
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
        if constexpr (fmm_order == 2) continue;

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
          if constexpr (fmm_order == 3) continue;

          // Octopole 
          res(m,q) += gc*H(i,q,m)*r[i]/d7;
          for(int j = 0 ; j < gdimension; ++j){
            // Octopole
            res(m,q) += -gc*3.5*(H(i,j,m)*r[q]+H(i,j,q)*r[m])*r[i]*r[j]/d9;
            #ifdef HEXA
            // Hexadecapole
            res(m,q) += gc*.5*X(i,j,q,m)*r[i]*r[j]/d9; 
            #endif 
            for(int k = 0 ; k < gdimension; ++k){
              // Octopole 
              res(m,q) += gc*7./6.*(9.*r[m]*r[q]/d2-(q==m))*
                H(i,j,k)*r[i]*r[j]*r[k]/d9; 
              #ifdef HEXA
              // Hexadecapole 
              res(m,q) += -gc*9./6.*(X(i,j,k,m)*r[q]+
                X(i,j,k,q)*r[m])*r[i]*r[j]*r[k]/d11;
              for(int l = 0; l < gdimension; ++l){
                // Hexadecapole 
                res(m,q) += gc*9./24.*(11.*r[m]*r[q]/d2-(q==m))*
                  X(i,j,k,l)*r[i]*r[j]*r[k]*r[l]/d11;
              }
              #endif // HEXA
            } // for k
          } // for j
        } // for i
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
    const point_t& dist_coordinates,
    const double& M, 
    const sym_tensor_rank4& X, 
    const sym_tensor_rank3& H, 
    const sym_tensor_rank2& Q)
  {

    if constexpr (fmm_order <= 2) return;

    if constexpr (fmm_order == 3) {
      double d = flecsi::distance(local_coordinates,dist_coordinates);
      point_t r = local_coordinates-dist_coordinates; 
      const double d2 = d*d; 
      const double d3 = d2*d; 
      const double d5 = d3*d2; 
      for(int m = 0; m < gdimension; ++m)
      for(int q = m; q < gdimension; ++q)
      for(int s = q; s < gdimension; ++s) {
        res(m,q,s) += gc*3.*M/d5*(
          (m==q)*r[s]+(q==s)*r[m]+(m==s)*r[q]-5.*r[m]*r[q]*r[s]/d2);
      }
      return;
    }

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
  }


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
    const sym_tensor_rank2& dfcdr = source->dfcdr();
    const sym_tensor_rank3& dfcdrdr = source->dfcdrdr();
    point_t cofm_coordinates = source->coordinates();

    point_t part_coordinates = sink->coordinates();
    point_t r = part_coordinates - cofm_coordinates;
    point_t grav = fc;
    double pot = pc;

    for(int i = 0 ; i < gdimension; ++i){
      pot += -r[i]*fc[i]; 
    }

    if constexpr (fmm_order > 1) {
      // The Jacobi
      for(int m=0;m<gdimension;++m){
        for(int i=0;i<gdimension;++i){
          grav[m] += dfcdr(m,i)*r[i];
          pot += -.5*r[m]*r[i]*dfcdr(m,i);
        } // for
      } // for

      if constexpr (fmm_order > 2) {
        // The hessian
        for(int m = 0; m < gdimension; ++m){
          for(int i = 0; i < gdimension; ++i){
            for(int j = 0 ; j < gdimension; ++j){
              grav[m] += .5*r[i]*r[j]*dfcdrdr(m,i,j); 
              pot += -1./6.*r[m]*r[i]*r[j]*dfcdrdr(m,i,j);
            } // for
          } // for
        } // for
      }
    }
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
        fmm_comms n = node_sources[k];
        gravitation_fc(pc, acc, p->coordinates(), n.coords, n.T,n.X,n.H,n.Q); 
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
    gravitation_fc(sink->pc(), sink->fc(), sink->coordinates(), 
      source->coordinates(), source->mass(), 
      source->hexa(), source->octo(), source->quad());
    if constexpr (fmm_order == 1) return;
    gravitation_dfcdr(sink->dfcdr(), sink->coordinates(), 
      source->coordinates(), source->mass(), 
      source->hexa(), source->octo(), source->quad());
    if constexpr (fmm_order == 2) return;
    gravitation_dfcdrdr(sink->dfcdrdr(), sink->coordinates(), 
      source->coordinates(), source->mass(), 
      source->hexa(), source->octo(), source->quad());
  }

  /**
  * @brief node<-particle interaction: update Taylor expansion coefficients
  */
  void taylor_p2c(
    node* sink,
    const type_t point_source_mass,
    const point_t& point_source_coords)
  {
    gravitation_fc(sink->pc(), sink->fc(), sink->coordinates(), 
      point_source_coords, point_source_mass, 0, 0, 0); 
    if constexpr (fmm_order == 1) return;
    gravitation_dfcdr(sink->dfcdr(), sink->coordinates(), 
      point_source_coords, point_source_mass, 0, 0, 0); 
    if constexpr (fmm_order == 2) return;
    gravitation_dfcdrdr(sink->dfcdrdr(), sink->coordinates(), 
      point_source_coords, point_source_mass, 0, 0, 0); 
  }

} // namespace fmm


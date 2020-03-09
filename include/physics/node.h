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
 * @file cofm.h
 * @author Julien Loiseau
 * @date January 2020
 * @brief Representation of a center of mass
 */

#pragma once

#include "space_vector.h"
#include "tree_topology/tree_types.h"
#include "user.h"
#include "tensor.h"

namespace flecsi {
  using sym_tensor_rank2 = flecsi::tensor_u<type_t,
                          symmetry_type::symmetric,gdimension,gdimension>;
  using sym_tensor_rank3 = flecsi::tensor_u<type_t, 
               symmetry_type::symmetric,gdimension,gdimension,gdimension>;
  using sym_tensor_rank4 = flecsi::tensor_u<type_t, 
    symmetry_type::symmetric,gdimension,gdimension,gdimension,gdimension>;
}

template<class KEY>
class node_u : public flecsi::topology::cofm_u<gdimension,type_t,KEY> {

  static const size_t dimension = gdimension;
  using element_t = type_t;
  using point_t = flecsi::space_vector_u<element_t, dimension>;
  using key_t = KEY; 

  using sym_tensor_rank2 = flecsi::sym_tensor_rank2;
  using sym_tensor_rank3 = flecsi::sym_tensor_rank3;
  using sym_tensor_rank4 = flecsi::sym_tensor_rank4;

public:

  node_u(): flecsi::topology::cofm_u<gdimension,type_t,KEY>()
  {
    X_ = {0}; H_ = {0}; Q_ = {0}; 
    pc_ = 0;
    fc_ = 0;
    dfcdr_ = 0;
    dfcdrdr_ = 0;
    affected_ = false;
  }

  node_u(const key_t & key): 
    flecsi::topology::cofm_u<gdimension,type_t,KEY>(key)
  {
    X_ = {0}; H_ = {0}; Q_ = {0}; 
    pc_ = 0;
    fc_ = 0;
    dfcdr_ = 0;
    dfcdrdr_ = 0;
    affected_ = false;
  }

  explicit node_u(const node_u& c): 
    flecsi::topology::cofm_u<gdimension,type_t,KEY>(c)
  {
    X_ = c.hexa(); H_ = c.octo(); Q_ = c.quad(); 
    pc_ = c.pc_; fc_ = c.fc_; dfcdr_ = c.dfcdr_; dfcdrdr_ = c.dfcdrdr_;
    affected_ = c.affected_;
  }

  const sym_tensor_rank4& hexa() const {return X_;}
  const sym_tensor_rank3& octo() const {return H_;}
  const sym_tensor_rank2& quad() const {return Q_;}

  sym_tensor_rank4& hexa() {return X_;}
  sym_tensor_rank3& octo() {return H_;}
  sym_tensor_rank2& quad() {return Q_;}

  const type_t& pc() const {return pc_;}
  const point_t& fc() const {return fc_;}
  const sym_tensor_rank2& dfcdr() const {return dfcdr_;}
  const sym_tensor_rank3& dfcdrdr() const {return dfcdrdr_;}

  type_t& pc() {return pc_;}
  point_t& fc() {return fc_;}
  sym_tensor_rank2& dfcdr() {return dfcdr_;}
  sym_tensor_rank3& dfcdrdr() {return dfcdrdr_;}

  void set_affected(const bool& affected) {affected_ = affected;} 
  bool affected() const {return affected_;}

private:
  sym_tensor_rank4 X_;
  sym_tensor_rank3 H_;
  sym_tensor_rank2 Q_; 

  type_t pc_;
  point_t fc_;
  sym_tensor_rank2 dfcdr_;
  sym_tensor_rank3 dfcdrdr_;

  bool affected_;

}; // class node


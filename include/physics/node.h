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

#ifndef node_h
#define node_h


#include "space_vector.h"
#include "tree_topology/tree_types.h"
#include "user.h"
#include "tensor.h"

template<class KEY>
class node_u : public flecsi::topology::cofm_u<gdimension,type_t,KEY> {

  static const size_t dimension = gdimension;
  using element_t = type_t;
  using point_t = flecsi::space_vector_u<element_t, dimension>;
  using key_t = KEY; 

public:

  node_u(): flecsi::topology::cofm_u<gdimension,type_t,KEY>(){}

  node_u(const key_t & key): 
    flecsi::topology::cofm_u<gdimension,type_t,KEY>(key){} 

  node_u(const node_u& c): 
    flecsi::topology::cofm_u<gdimension,type_t,KEY>(c)
  {
    //X_ = c.hexa(); 
    //H_ = c.octo(); 
    //Q_ = c.quad(); 
  }

  //const flecsi::tensor_u<double, symmetry_type::symmetric,3,3,3,3>& hexa() const
  //{return X_;}
  //const flecsi::tensor_u<double, symmetry_type::symmetric,3,3,3>& octo() const
  //{return H_;}
  //const flecsi::tensor_u<double, symmetry_type::symmetric,3,3>& quad() const
  //{return Q_;}

  flecsi::tensor_u<double, symmetry_type::symmetric,3,3,3,3>& hexa()
  {return X_;}
  flecsi::tensor_u<double, symmetry_type::symmetric,3,3,3>& octo()
  {return H_;}
  flecsi::tensor_u<double, symmetry_type::symmetric,3,3>& quad()
  {return Q_;}

private:
  flecsi::tensor_u<double, symmetry_type::symmetric,3,3,3,3> X_;
  flecsi::tensor_u<double, symmetry_type::symmetric,3,3,3> H_;
  flecsi::tensor_u<double, symmetry_type::symmetric,3,3> Q_; 

}; // class node

#endif // node_h

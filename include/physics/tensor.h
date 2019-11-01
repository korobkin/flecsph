/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2019 Triad National Security, LLC
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
 * @file tensor.h
 * @author Oleg Korobkin
 * @date November 2019
 * @brief Tensors rank2 class
 */

#ifndef TENSOR_H
#define TENSOR_H

#include <flecsi/geometry/point.h>
#include <flecsi/utils/common.h>
#include <flecsi/utils/dimensioned_array.h>

//----------------------------------------------------------------------------//
//! Enumeration for axes.
//----------------------------------------------------------------------------//

enum class axis : size_t { x = 0, y = 1, z = 2 };

namespace tensor_indices {
  enum dimensional_axes { xx = 0, xy = 1, xz = 2,
                          yx = 3, yy = 4, yz = 5,
                          zx = 6, zy = 7, zz = 8};
}

namespace flecsi {

/*!
  \class generic_tensor_rank2 tensor.h
  \brief This class defines an interface for operations on generic tensors of
  rank 2 in D-dimensional space. Does not assume any symmetries.

  The generic_tensor_rank2 type is implemented using \ref dimensioned_array.
  Look there for more information on the vector_t interface.

  \tparam T   Data type
  \tparam D   Dimension of the vector space on which tensor is defined
  \tparam N   The namespace of the array.  This is a dummy parameter
              that is useful for creating distinct types that alias
              dimensioned_array_u.
 */
template<typename T, size_t D, size_t N>
struct generic_tensor_rank2_u
  : public utils::dimensioned_array_u<T, D*D, N> {

  using parent_class = utils::dimensioned_array_u<T, D*D, N>;

  /*!
   \function  operator[](const size_t & ind)
   \brief     Data write operator with enumerated type access:
              A[xy] or B[zz], etc. See 'tensor_indices' enum above

   \tparam    T    data type
   \param[in] ind  flattened index, from 0 to D^2-1
   \return    reference to data element at the index location
   */
  T & operator[](const size_t & ind) {
    return parent_class::operator[](ind);
  } // operator []


  /*!
   \function  operator()(const size_t & ind)
   \brief     Data read operator (accessor) via single-index
              enumerated access: A[xy] or B[zz], etc.
              See 'tensor_indices' enum above

   \tparam    T    data type
   \param[in] ind  flattened index, from 0 to D^2-1
   \return    data value of the element at the index location
   */
  T const operator()(const size_t & ind) const {
    return parent_class::operator[](ind);
  } // operator ()


  /*!
   \function  operator()(const size_t & i, const size_t & j)
   \brief     Data write operator (mutator) via two-index notation:
              A(1,2) = 3, x = B(3,3) etc.

   \tparam    T     data type
   \param[in] i     first index, from 0 to D-1
   \param[in] j     second index, from 0 to D-1
   \return    data  reference to the element at the index location
   */
  T & operator()(const size_t & i, const size_t & j) {
    return operator[](i*D + j);
  } // operator ()


  /*!
   \function  operator()(const size_t & i, const size_t & j)
   \brief     Data read operator (accessor) via two-index notation:
              q_xy = Q(1,2) etc. Suitable for accessing elements of
              a tensor which was defined as a const.

   \tparam    T     data type
   \param[in] i     first index, from 0 to D-1
   \param[in] j     second index, from 0 to D-1
   \return    data  reference to the element at the index location
   */
  T const operator()(const size_t & i, const size_t & j) const {
    return operator()(i*D + j);
  } // operator ()

}; // struct generic_tensor_rank2


/*!
  \function      operator<<(std::ostream, generic_tensor_rank2)
  \brief         Output stream operator for the generic tensor

  \tparam T      Data type
  \tparam D      Dimension of the vector space on which tensor is defined
  \tparam N      FleCSI namespace

  \param stream  The output stream.
  \param a       The tensor to output
 */
template<typename T, size_t D, size_t N>
std::ostream &
operator<<(std::ostream & stream,
  generic_tensor_rank2_u<T, D, N> const & a) {
  stream << "[";
  for(size_t i = 0; i < D; i++) {
    stream << "[" << a(i, 0);
    for(size_t j = 1; j < D; j++) {
      stream << ", " << a(i,j);
    }
    stream << "]";
    if (i + 1 < D) stream << ", ";
  } // for
  stream << " ]";

  return stream;
} // operator << generic_tensor_rank2


/*!
  \class symmetric_tensor_rank2 tensor.h
  \brief Interface for symmetric tensors of rank 2 in D-dimensional space.
         Inherits basic arithmetic operations from \ref dimensioned_array.

  \tparam T   Data type
  \tparam D   Dimension of the vector space on which tensor is defined
  \tparam N   FleCSI namespace
 */
template<typename T, size_t D, size_t N>
struct symmetric_tensor_rank2_u
: public utils::dimensioned_array_u<T, D*(D+1)/2, N> {

  using parent_class = utils::dimensioned_array_u<T, D*(D+1)/2, N>;

  // tensor size: diagonal + upper (or lower) triangle
  static constexpr size_t TENSOR_SIZE = D*(D + 1)/2;

  // index mappings for 2D and 3D
  static constexpr size_t map_2d[4] = { 0, 1,
                                        1, 2 };

  static constexpr size_t map_3d[9] = { 0, 1, 2,
                                        1, 3, 4,
                                        2, 4, 5 };

  /*!
   \function  operator[](const size_t & ind)
   \brief     Data write operator with enumerated type access:
              A[xy] or B[zz], etc. All index pairs are ok.

   \tparam    T    data type
   \param[in] ind  flattened index, from 0 to D^2-1
   \return    reference to data element at the index location
   */
  T & operator[](const size_t & ind) {
    if constexpr (D == 1) {
      assert (ind == 0);
      return parent_class::operator[](0);
    }

    int k;
    assert (ind >=0 && ind < D*D);
    if constexpr (D == 2)
      k = map_2d[ind];

    if constexpr (D == 3)
      k = map_3d[ind];

    return parent_class::operator[](k);
  } // operator []


  /*!
   \function  operator()(const size_t & ind)
   \brief     Data read operator (accessor) via single-index
              enumerated access: A[xy] or B[zz], etc.

   \tparam    T    data type
   \param[in] ind  flattened index, from 0 to D^2-1
   \return    data value of the element at the index location
   */
  T operator()(const size_t & ind) const {
    if constexpr (D == 1) {
      assert (ind == 0);
      return parent_class::operator[](0);
    }

    int k = -1;
    assert (ind >=0 && ind < D*D);
    if constexpr (D == 2)
      k = map_2d[ind];

    if constexpr (D == 3)
      k = map_3d[ind];

    return parent_class::operator[](k);
  } // operator ()


  /*!
   \function  operator()(const size_t & i, const size_t & j)
   \brief     Data write operator (mutator) via two-index notation:
              A(1,2) = 3, x = B(3,3) etc.

   \tparam    T     data type
   \param[in] i     first index, from 0 to D-1
   \param[in] j     second index, from 0 to D-1
   \return    data  reference to the element at the index location
   */
  T & operator()(const size_t & i, const size_t & j) {
    return operator[](i*D + j);
  } // operator ()


  /*!
   \function  operator()(const size_t & i, const size_t & j)
   \brief     Data read operator (accessor) via two-index notation:
              q_xy = Q(1,2) etc. Suitable for accessing elements of
              a tensor which was defined as a const.

   \tparam    T     data type
   \param[in] i     first index, from 0 to D-1
   \param[in] j     second index, from 0 to D-1
   \return    data  reference to the element at the index location
   */
  T const operator()(const size_t & i, const size_t & j) const {
    return operator()(i*D + j);
  } // operator ()

}; // struct symmetric_tensor_rank2

/*!
  \function      operator<<(std::ostream, symmetric_tensor_rank2)
  \brief         Output stream operator for the symmetric tensor

  \tparam T      Data type
  \tparam D      Dimension of the vector space on which tensor is defined
  \tparam N      FleCSI namespace

  \param stream  The output stream.
  \param a       The symmetric tensor to output
 */
template<typename T, size_t D, size_t N>
std::ostream &
operator<<(std::ostream & stream,
  symmetric_tensor_rank2_u<T, D, N> const & a) {
  stream << "[";
  for(size_t i = 0; i < D; i++) {
    stream << "[" << a(i, 0);
    for(size_t j = 1; j < D; j++) {
      stream << ", " << a(i,j);
    }
    stream << "]";
    if (i + 1 < D) stream << ", ";
  } // for
  stream << " ]";

  return stream;
} // operator << symmetric_tensor_rank2


/*!
  \class antisymmetric_tensor_rank2 tensor.h
  \brief Interface for anti-symmetric rank-2 tensors in D-dimensional space.
         Inherits basic arithmetic operations from \ref dimensioned_array.

  \tparam T   Data type
  \tparam D   Dimension of the vector space on which tensor is defined
  \tparam N   FleCSI namespace
 */
template<typename T, size_t D, size_t N>
struct antisymmetric_tensor_rank2_u
: public utils::dimensioned_array_u<T, D*(D-1)/2, N> {

  using parent_class = utils::dimensioned_array_u<T, D*(D-1)/2, N>;

  // tensor size: upper triangle of the tensor matrix
  static const size_t TENSOR_SIZE = D*(D - 1)/2;

  // index mappings for 2D and 3D
  static constexpr int map_2d[4] = {  0,  1,
                                     -1,  0};

  static constexpr int map_3d[9] = {  0, 1, 2,
                                     -1, 0, 3,
                                     -2,-3, 0};

  /*!
   \function  operator[](const size_t & ind)
   \brief     Data write operator with enumerated type access:
              A[xy] or B[xz], etc.
              Valid indices:
              - in 1D: none (no storage in 1D);
              - 2D; xy only;
              - 3D: xy, xz, yz.
              If invalid index is given, assertion fault is produced.

   \tparam    T    data type
   \param[in] ind  flattened index: see acceptable values above
   \return    reference to data element at the index location.
              Notice how lower-triangle and diagonal elements are
              not accessible: the former are opposite of upper-triangle,
              while the latter are zeros.
   */
  T & operator[](const size_t & ind) {
    int k = -1;
    if constexpr (D == 2)
       k = map_2d[ind];

    if constexpr (D == 3)
       k = map_3d[ind];

    assert (k > 0);
    return parent_class::operator[](k-1);

  } // operator []


  /*!
   \function  operator()(const size_t & ind)
   \brief     Data read operator (accessor) via single-index
              enumerated access: A[xy] or B[zz], etc.
              For data reading, all indices (xx, xy, .. zz)
              are OK.

   \tparam    T    data type
   \param[in] ind  flattened index, from 0 to D^2-1
   \return    data value of the element at the index location
   */
  T const operator()(const size_t & ind) const {
    int k = -1;
    if constexpr (D == 2)
       k = map_2d[ind];

    if constexpr (D == 3)
       k = map_3d[ind];

    if (k == 0) return (T)0;
    int sgn = (k>0) - (k<0);
    k = abs(k) - 1;
    assert (k >= 0 and k < TENSOR_SIZE);
    return sgn*(parent_class::operator[](k));
  } // operator ()


  /*!
   \function   operator()(const size_t & i, const size_t & j)
   \brief      Data read operator (accessor) via two-index notation:
               q_xy = Q(1,2) etc. Suitable for accessing elements of
               a tensor which was defined as a const.
               Notice that a two-index mutator for antisymmetric
               tensor is not defined. Use a single-index mutator 
               instead.

   \tparam     T     data type
   \param[in]  i     first index, from 0 to D-1
   \param[in]  j     second index, from 0 to D-1
   \return     data  reference to the element at the index location
   */
  T const operator()(const size_t & i, const size_t & j) const {
    return operator()(i*D + j);
  } // operator ()

}; // struct antisymmetric_tensor_rank2


/*!
  \function      operator<<(std::ostream, antisymmetric_tensor_rank2)
  \brief         Output stream operator for the antisymmetric tensor

  \tparam T      Data type
  \tparam D      Dimension of the vector space on which tensor is defined
  \tparam N      FleCSI namespace

  \param stream  The output stream.
  \param a       The antisymmetric tensor to output
 */
template<typename T, size_t D, size_t N>
std::ostream &
operator<<(std::ostream & stream,
  antisymmetric_tensor_rank2_u<T, D, N> const & a) {
  stream << "[";
  for(size_t i = 0; i < D; i++) {
    stream << "[" << a(i, 0);
    for(size_t j = 1; j < D; j++) {
      stream << ", " << a(i,j);
    }
    stream << "]";
    if (i + 1 < D) stream << ", ";
  } // for
  stream << " ]";

  return stream;
} // operator << antisymmetric_tensor_rank2

} // namespace flecsi

#endif // TENSOR_H

/*
// Usage example
#include <iostream>

int main() {

  using namespace std;
  using namespace tensor_indices;

  cout << "--- Antisymmetric tensor: ---" << endl;
  using antitensor_t = antisymmetric_tensor_rank2<double, 3, 0>;
  antitensor_t A{0};
  A[xy] = 1.2; A[xy] = 2.3; A[yz] = 1.3;
  cout << "size of A: " << A.size() << endl;
  cout << "tensor A: "  << endl << A << endl;

  cout << "--- Symmetric tensor: ---" << endl;
  using sym_tensor_t = symmetric_tensor_rank2<double, 3, 0>;
  sym_tensor_t S{0};
  S[xy] = 1.2; S[yz] = 2.3; S[xz] = 1.3;
  S(1,1) = 2.0; S[xx] = 0.1;
  cout << "size of S: " << S.size() << endl;
  cout << "tensor S: "  << endl << S << endl;

  cout << "--- Generic tensor: ---" << endl;
  using gen_tensor_t = generic_tensor_rank2<double, 3, 0>;
  gen_tensor_t G{0};
  G[xy] = 1.2; G[yz] = 2.3; G[xz] = 1.3;
  G(1,1) = 2.0; G[xx] = 0.1;
  cout << "size of G: " << G.size() << endl;
  cout << "tensor G: "  << endl << G << endl;

}
*/

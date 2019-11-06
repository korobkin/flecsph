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
 * @brief General tensors class
 */

#ifndef TENSOR_H
#define TENSOR_H

#include <array>
#include <cmath>
#include <ostream>

/*!
 \class symmetry_type tensor.h
 \brief Defines symmetry type for rank-2 tensors
 */
enum class symmetry_type {
  generic,
  symmetric
};

//----------------------------------------------------------------------------//
//! Enumeration for axes.
//----------------------------------------------------------------------------//
namespace tensor_indices {
  enum rank1_index { _x = 0, _y = 1, _z = 2, _t = 3 };
  enum rank2_index { xx = 0, xy = 1, xz = 2, xt = 3,
                     yx = 4, yy = 5, yz = 6, yt = 7,
                     zx = 8, zy = 9, zz =10, zt =11,
                     tx =12, ty =13, tz =14, tt =15 };

} // namespace tensor_indices

namespace flecsph {

/*!
  \class tensor_u tensor.h
  \brief This class defines an interface for operations on generic tensors of
         arbitrary rank in an arbitrary product of vector spaces over the field
         of type T.

  \tparam T      Data field: e.g. float, double, complex...
  \tparam ST     Symmetry type. So far, only generic and symmetric for rank-2 
                 tensors have been implemented
  \tparam Ds..   Dimensions of the product vector space where the tensor is 
                 acting
 */
template <class T, symmetry_type ST, auto... Ds>   // variadic recursion base
struct tensor_u {
  static constexpr auto size() { return 1; };      // data size
  static constexpr auto RANK = 0;                  // tensor rank
  static constexpr size_t DIM[0] = {};             // tensor dimensions
  static constexpr auto multiindex() { return 0; }
};

template <class T, symmetry_type ST, auto D, auto... Ds>
struct tensor_u<T, ST, D, Ds...> {

  // data size
  static constexpr auto size() {
    if constexpr (ST == symmetry_type::generic)
      return D * tensor_u<T, ST, Ds...>::size();

    if constexpr (ST == symmetry_type::symmetric) // rank 2 only
      return D*(D+1)/2;

  }

  // tensor rank
  static constexpr auto RANK = 1 + tensor_u<T, ST, Ds...>::RANK;

  // tensor dimensions
  static constexpr size_t DIM[RANK] = {D, Ds...};

  // access flattened data via integer-type index
  T& operator[](const size_t ind) { 
    return data_[ind]; 
  }

  // for tensors of rank 1: access using A[_x] etc. notation
  T& operator[](const enum tensor_indices::rank1_index ind) {
    assert (RANK == 1);
    assert (ind >= 0 and ind < size());
    return data_[ind];
  }

  // tensors of rank 2: access using A[yx], A[zz] etc. notation
  T& operator[](const enum tensor_indices::rank2_index ind) {
    assert (RANK == 2);
    using namespace tensor_indices;

    // generic case: no symmetry
    if constexpr (ST == symmetry_type::generic) {
      if constexpr (D == 1) {
        assert (ind == 0);
        return data_[0];
      }
      if constexpr (D == 2) {
        constexpr int remap[] = { 0, 2, -1, -1,
                                  1, 3};
        assert (ind < 6);
        assert (remap[ind] >= 0);
        return data_[remap[ind]];
      }
      if constexpr (D == 3) {
        constexpr int remap[] = { 0, 3, 6, -1,
                                  1, 4, 7, -1,
                                  2, 5, 8 };
        assert (ind < 11);
        assert (remap[ind] >= 0);
        return data_[remap[ind]];
      }
      if constexpr (D == 4) {
        constexpr int remap[] = { 0, 4,  8, 12,
                                  1, 5,  9, 13,
                                  2, 6, 10, 14,
                                  3, 7, 11, 15};
        return data_[remap[ind]];
      }
    }

    // symmetric tensor of rank 2
    if constexpr (ST == symmetry_type::symmetric) {
      if constexpr (D == 1) {
        assert (ind == 0);
        return data_[0];
      }
      if constexpr (D == 2) {
        constexpr int remap[] = { 0, 1, -1, -1,
                                  1, 2};
        assert (ind < 6);
        assert (remap[ind] >= 0);
        return data_[remap[ind]];
      }
      if constexpr (D == 3) {
        constexpr int remap[] = { 0, 1, 3, -1,
                                  1, 2, 4, -1,
                                  3, 4, 5 };
        assert (ind < 11);
        assert (remap[ind] >= 0);
        return data_[remap[ind]];
      }
      if constexpr (D == 4) {
        constexpr int remap[] = { 0, 1, 3, 6,
                                  1, 2, 4, 7,
                                  3, 4, 5, 8,
                                  6, 7, 8, 9};

        return data_[remap[ind]];
      }
    }

    assert (false);
  }

  // multi-index function: converts multi-index into flat index
  // constexpr for compile-time eval
  template <class Ind, class... Inds>
  static constexpr auto multiindex(Ind&& i, Inds&&... inds) {
    if constexpr (ST == symmetry_type::generic)
      return D * tensor_u<T, ST, Ds...>::multiindex(inds...) + i;

    if constexpr (ST == symmetry_type::symmetric)
      return (((inds+...) + 1)*(inds+...)/2 + i) > ((inds+...) + i*(i+1)/2) ?
             (((inds+...) + 1)*(inds+...)/2 + i) : ((inds+...) + i*(i+1)/2);

  }

  // multi-index access interface
  // - constexpr:    for compile-time eval
  // decltype(auto): for latest-time type evaluation (e.g. value or ref
  //                 depending on call context)
  template <class... Inds>
  constexpr decltype(auto) operator()(Inds&&... inds) {
    return data_[multiindex(inds...)];
  }


private:
  T data_[size()];

}; // tensor_u


template<typename T, size_t D>
using space_vector_u = tensor_u<T, symmetry_type::generic, D>;

} // namespace flecsph

#endif // TENSOR_H


#if 0
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

#endif // DEBUG

/*
// Usage example
#include <iostream>

int main() {

  using namespace std;
  using namespace flecsph;
  using namespace tensor_indices;

  cout << "--- Vectors: -------------------" << endl;
  using vector_1d_t = space_vector<double, 1>;
  vector_1d_t b1;
  cout << "1D: ";
  b1[0]  = 1.8; cout << "b1[0]  = " << b1[0]  << "\n";
  cout << "    ";
  b1[_x] = 1.9; cout << "b1[_x] = " << b1[_x] << "\n";
  cout << "    ";
  b1(0)  = 2.1; cout << "b1(0)  = " << b1(0)  << "\n";

  using vector_2d_t = space_vector<double, 2>;
  vector_2d_t b2;
  cout << "2D: ";
  b2[0]  = 1.8; b2[1]  = -1.8;
  cout << "{b2[0],  b2[1] } = {" << b2[0]  << ", " << b2[1] << "}\n";
  cout << "    ";
  b2[_x] = 1.9; b2[_y] = -1.9;
  cout << "{b2[_x], b2[_y]} = {" << b2[_x] << ", " << b2[_y] << "}\n";
  cout << "    ";
  b2(0)  = 2.1; b2(1)  = -2.1;
  cout << "{b2(0),  b2(1) } = {" << b2(0)  << ", " << b2(1) << "}\n";

  using vector_3d_t = space_vector<double, 3>;
  vector_3d_t b3;
  cout << "3D: ";
  b3[0]  = 1.1; b3[1]  = 1.2; b3[2] = 1.3;
  cout <<"{b3[0 .. 2]} = {"<<b3[0]<<", "<<b3[1]<<", "<<b3[2]<< "}\n";
  cout << "    ";
  b3[_x] = 2.1; b3[_y] = 2.2; b3[_z] = 2.3;
  cout << "{b3[_x.._z]} = {"<<b3[_x]<<", "<<b3[_y]<<", "<<b3[_z]<< "}\n";
  cout << "    ";
  b3(0)  = 3.1; b3(1)  = 3.2; b3(2)  = 3.3;
  cout <<"{b3(0 .. 2)} = {"<<b3(0)<<", "<<b3(1)<<", "<<b3(2)<< "}\n";

  //cout << "--- Antisymmetric tensor: ---" << endl;
  //using antitensor_t = antisymmetric_tensor_rank2<double, 3, 0>;
  //antitensor_t A{0};
  //A[xy] = 1.2; A[xy] = 2.3; A[yz] = 1.3;
  //cout << "size of A: " << A.size() << endl;
  //cout << "tensor A: "  << endl << A << endl;

  //cout << "--- Symmetric tensor: ---" << endl;
  //using sym_tensor_t = symmetric_tensor_rank2<double, 3, 0>;
  //sym_tensor_t S{0};
  //S[xy] = 1.2; S[yz] = 2.3; S[xz] = 1.3;
  //S(1,1) = 2.0; S[xx] = 0.1;
  //cout << "size of S: " << S.size() << endl;
  //cout << "tensor S: "  << endl << S << endl;

  //cout << "--- Generic tensor: ---" << endl;
  //using gen_tensor_t = generic_tensor_rank2<double, 3, 0>;
  //gen_tensor_t G{0};
  //G[xy] = 1.2; G[yz] = 2.3; G[xz] = 1.3;
  //G(1,1) = 2.0; G[xx] = 0.1;
  //cout << "size of G: " << G.size() << endl;
  //cout << "tensor G: "  << endl << G << endl;

}
*/

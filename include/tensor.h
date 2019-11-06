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
#include <flecsi/utils/common.h>

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

namespace flecsi {

template<typename... CONDITIONS>
struct and_ : std::true_type {};

template<typename CONDITION, typename... CONDITIONS>
struct and_<CONDITION, CONDITIONS...>
  : std::conditional<CONDITION::value, and_<CONDITIONS...>, std::false_type>::
      type {}; // struct and_

template<typename TARGET, typename... TARGETS>
using are_type_u = and_<std::is_same<TARGETS, TARGET>...>;

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

  //! Default constructor.
  tensor_u() = default;

  //! Default copy constructor.
  tensor_u(tensor_u const &) = default;

  //! Constructor (fill with given value).
  tensor_u(T const & val) {
    for (size_t i=0; i<size(); ++i)
      data_[i]= val;
  }

  //! Initializer list constructor.
  tensor_u(std::initializer_list<T> list) {
    assert(list.size() == size() && "dimension size mismatch");
    size_t i = 0;
    for(auto it = list.begin(); it!=list.end(); ++it, ++i)
      data_[i]= *it;
  }

  //! Variadic constructor.
  template<class... ARGS>
  tensor_u(T arg, ARGS... args) {
    std::initializer_list<T> list = {arg, args...};
    assert(list.size() == size() && "dimension size mismatch");
    size_t i = 0;
    for(auto it = list.begin(); it!=list.end(); ++it, ++i)
      data_[i]= *it;
  } // tensor_u


  //! Assignment operator.
  tensor_u& operator=(tensor_u const & rhs) {
    if(this != &rhs) {
      for (size_t i=0; i<size(); ++i)
        data_[i]= rhs.data_[i];
    } // if

    return *this;
  } // operator =

  //! Assignment operator with a scalar rhs
  tensor_u& operator=(const T & val) {
    for (size_t i=0; i<size(); ++i)
      data_[i]= val;
    return *this;
  } // operator =

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
  constexpr decltype(auto) operator[](const size_t ind) {
    return data_[ind];
  }

  // access data in const instance
  constexpr decltype(auto) operator[](const size_t ind) const {
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

  //--------------------------------------------------------------------------//
  // Macro to avoid code replication.
  //--------------------------------------------------------------------------//

#define define_operator(op)                                                    \
  tensor_u & operator op(tensor_u const & rhs) {                               \
    if(this != &rhs) {                                                         \
      for(size_t i{0}; i < size(); i++) {                                      \
        data_[i] op rhs[i];                                                    \
      } /* for */                                                              \
    } /* if */                                                                 \
                                                                               \
    return *this;                                                              \
  }

  //--------------------------------------------------------------------------//
  // Macro to avoid code replication.
  //--------------------------------------------------------------------------//

#define define_operator_type(op)                                               \
  tensor_u & operator op(T val) {                                              \
    for(size_t i{0}; i < size(); i++) {                                        \
      data_[i] op val;                                                         \
    } /* for */                                                                \
                                                                               \
    return *this;                                                              \
  }

  //--------------------------------------------------------------------------//
  //! Addition/Assignment operator.
  //--------------------------------------------------------------------------//

  define_operator(+=);

  //--------------------------------------------------------------------------//
  //! Addition/Assignment operator.
  //--------------------------------------------------------------------------//

  define_operator_type(+=);

  //--------------------------------------------------------------------------//
  //! Subtraction/Assignment operator.
  //--------------------------------------------------------------------------//

  define_operator(-=);

  //--------------------------------------------------------------------------//
  //! Subtraction/Assignment operator.
  //--------------------------------------------------------------------------//

  define_operator_type(-=);

  //--------------------------------------------------------------------------//
  //! Division/Assignment operator.
  //--------------------------------------------------------------------------//

  define_operator_type(*=);

  //--------------------------------------------------------------------------//
  //! Division/Assignment operator.
  //--------------------------------------------------------------------------//

  define_operator_type(/=);

  //! \brief Division operator involving a constant.
  //! \param[in] val The constant on the right hand side of the operator.
  //! \return A reference to the current object.
  tensor_u operator/(T val) {
    tensor_u tmp(*this);
    tmp /= val;

    return tmp;
  } // operator /


private:
  T data_[size()];

}; // tensor_u


/*!
  \function      operator+(std::ostream, tensor_u)
  \brief         Addition operator between two tensors

  \tparam T      Data type
  \tparam ST     Symmetry class
  \tparam Ds...  Dimensions of the vector space on which tensor is defined

  \param a       first tensor
  \param b       second tensor
 */
template <class T, symmetry_type ST, auto... Ds>
tensor_u<T, ST, Ds...>
operator+(
    const tensor_u<T, ST, Ds...> & a,
    const tensor_u<T, ST, Ds...> & b) {
  tensor_u<T, ST, Ds...> tmp(a);
  tmp += b;
  return tmp;
} // operator +


/*!
  \function      operator-(std::ostream, tensor_u)
  \brief         Subtraction operator between two tensors

  \tparam T      Data type
  \tparam ST     Symmetry class
  \tparam Ds...  Dimensions of the vector space on which tensor is defined

  \param a       first tensor
  \param b       second tensor
 */
template <class T, symmetry_type ST, auto... Ds>
tensor_u<T, ST, Ds...>
operator-(
    const tensor_u<T, ST, Ds...> & a,
    const tensor_u<T, ST, Ds...> & b) {
  tensor_u<T, ST, Ds...> tmp(a);
  tmp -= b;
  return tmp;
} // operator +


/*!
  \function      operator<<(std::ostream, tensor_u)
  \brief         Output stream operator for the generic tensor

  \tparam T      Data type
  \tparam ST     Symmetry class
  \tparam Ds...  Dimensions of the vector space on which tensor is defined

  \param stream  The output stream.
  \param a       The tensor to output
 */
template <class T, symmetry_type ST, auto... Ds>
std::ostream &
operator<<(std::ostream & stream,
    tensor_u<T, ST, Ds...> const & a) {
  using tensor = tensor_u<T, ST, Ds...>;
  if constexpr (tensor::RANK == 0) {
    stream << "[]";
  }
  else {
    stream << "[" << a[0];
    for(size_t i = 1; i < tensor::size(); ++i) {
      stream << ", " << a[i];
    }
    stream << "]";
  }
  return stream;
} // operator << tensor_u


/*!
  \function      operator*(tensor_u, T)
  \brief         Scalar multiplication operator for tensors

  \tparam T      Data type
  \tparam ST     Symmetry class
  \tparam Ds...  Dimensions of the vector space on which tensor is defined

  \param X       Tensor
  \param a       scalar
 */
template <class T, symmetry_type ST, auto... Ds>
tensor_u<T, ST, Ds...>
operator*(const tensor_u<T, ST, Ds...> & X, const T & a) {
  tensor_u<T, ST, Ds...> tmp(X);
  tmp *= a;
  return tmp;
} // operator *

template <class T, symmetry_type ST, auto... Ds>
tensor_u<T, ST, Ds...>
operator*(const T & a, const tensor_u<T, ST, Ds...> & X) {
  tensor_u<T, ST, Ds...> tmp(X);
  tmp *= a;
  return tmp;
} // operator *

/*!
  \function      operator*(tensor_u, T)
  \brief         Divide tensor by a scalar

  \tparam T      Data type
  \tparam ST     Symmetry class
  \tparam Ds...  Dimensions of the vector space on which tensor is defined

  \param X       Tensor
  \param a       scalar
 */
template <class T, symmetry_type ST, auto... Ds>
tensor_u<T, ST, Ds...>
operator/(const tensor_u<T, ST, Ds...> & X, const T & a) {
  tensor_u<T, ST, Ds...> tmp(X);
  tmp /= a;
  return tmp;
} // operator /

} // namespace flecsi

#endif // TENSOR_H


/*
// Usage example
#include <iostream>

int main() {

  using namespace std;
  using namespace flecsph;
  using namespace tensor_indices;


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

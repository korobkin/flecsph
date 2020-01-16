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
  enum rank3_index { xxx = 0, xxy = 1, xxz = 2, xxt = 3,
                     xyx = 4, xyy = 5, xyz = 6, xyt = 7,
                     xzx = 8, xzy = 9, xzz =10, xzt =11,
                     xtx =12, xty =13, xtz =14, xtt =15,
                     yxx =16, yxy =17, yxz =18, yxt =19,
                     yyx =20, yyy =21, yyz =22, yyt =23,
                     yzx =24, yzy =25, yzz =26, yzt =27,
                     ytx =28, yty =29, ytz =30, ytt =31,
                     zxx =32, zxy =33, zxz =34, zxt =35,
                     zyx =36, zyy =37, zyz =38, zyt =39,
                     zzx =40, zzy =41, zzz =42, zzt =43,
                     ztx =44, zty =45, ztz =46, ztt =47,
                     txx =48, txy =49, txz =50, txt =51,
                     tyx =52, tyy =53, tyz =54, tyt =55,
                     tzx =56, tzy =57, tzz =58, tzt =59,
                     ttx =60, tty =61, ttz =62, ttt =63 };
  enum rank4_index { xxxx =  0, xxxy =  1, xxxz =  2, xxxt =  3,
                     xxyx = 4, xxyy =  5, xxyz =  6, xxyt =  7,
                     xxzx = 8, xxzy =  9, xxzz = 10, xxzt = 11,
                     xxtx = 12, xxty = 13, xxtz = 14, xxtt = 15,
                     xyxx = 16, xyxy = 17, xyxz = 18, xyxt = 19,
                     xyyx = 20, xyyy = 21, xyyz = 22, xyyt = 23,
                     xyzx = 24, xyzy = 25, xyzz = 26, xyzt = 27,
                     xytx = 28, xyty = 29, xytz = 30, xytt = 31,
                     xzxx = 32, xzxy = 33, xzxz = 34, xzxt = 35,
                     xzyx = 36, xzyy = 37, xzyz = 38, xzyt = 39,
                     xzzx = 40, xzzy = 41, xzzz = 42, xzzt = 43,
                     xztx = 44, xzty = 45, xztz = 46, xztt = 47,
                     xtxx = 48, xtxy = 49, xtxz = 50, xtxt = 51,
                     xtyx = 52, xtyy = 53, xtyz = 54, xtyt = 55,
                     xtzx = 56, xtzy = 57, xtzz = 58, xtzt = 59,
                     xttx = 60, xtty = 61, xttz = 62, xttt = 63,
                     yxxx = 64, yxxy = 65, yxxz = 66, yxxt = 67,
                     yxyx = 68, yxyy = 69, yxyz = 70, yxyt = 71,
                     yxzx = 72, yxzy = 73, yxzz = 74, yxzt = 75,
                     yxtx = 76, yxty = 77, yxtz = 78, yxtt = 79,
                     yyxx = 80, yyxy = 81, yyxz = 82, yyxt = 83,
                     yyyx = 84, yyyy = 85, yyyz = 86, yyyt = 87,
                     yyzx = 88, yyzy = 89, yyzz = 90, yyzt = 91,
                     yytx = 92, yyty = 93, yytz = 94, yytt = 95,
                     yzxx = 96, yzxy = 97, yzxz = 98, yzxt = 99,
                     yzyx =100, yzyy =101, yzyz =102, yzyt =103,
                     yzzx =104, yzzy =105, yzzz =106, yzzt =107,
                     yztx =108, yzty =109, yztz =110, yztt =111,
                     ytxx =112, ytxy =113, ytxz =114, ytxt =115,
                     ytyx =116, ytyy =117, ytyz =118, ytyt =119,
                     ytzx =120, ytzy =121, ytzz =122, ytzt =123,
                     yttx =124, ytty =125, yttz =126, yttt =127,
                     zxxx =128, zxxy =129, zxxz =130, zxxt =131,
                     zxyx =132, zxyy =133, zxyz =134, zxyt =135,
                     zxzx =136, zxzy =137, zxzz =138, zxzt =139,
                     zxtx =140, zxty =141, zxtz =142, zxtt =143,
                     zyxx =144, zyxy =145, zyxz =146, zyxt =147,
                     zyyx =148, zyyy =149, zyyz =150, zyyt =151,
                     zyzx =152, zyzy =153, zyzz =154, zyzt =155,
                     zytx =156, zyty =157, zytz =158, zytt =159,
                     zzxx =160, zzxy =161, zzxz =162, zzxt =163,
                     zzyx =164, zzyy =165, zzyz =166, zzyt =167,
                     zzzx =168, zzzy =169, zzzz =170, zzzt =171,
                     zztx =172, zzty =173, zztz =174, zztt =175,
                     ztxx =176, ztxy =177, ztxz =178, ztxt =179,
                     ztyx =180, ztyy =181, ztyz =182, ztyt =183,
                     ztzx =184, ztzy =185, ztzz =186, ztzt =187,
                     zttx =188, ztty =189, zttz =190, zttt =191,
                     txxx =192, txxy =193, txxz =194, txxt =195,
                     txyx =196, txyy =197, txyz =198, txyt =199,
                     txzx =200, txzy =201, txzz =202, txzt =203,
                     txtx =204, txty =205, txtz =206, txtt =207,
                     tyxx =208, tyxy =209, tyxz =210, tyxt =211,
                     tyyx =212, tyyy =213, tyyz =214, tyyt =215,
                     tyzx =216, tyzy =217, tyzz =218, tyzt =219,
                     tytx =220, tyty =221, tytz =222, tytt =223,
                     tzxx =224, tzxy =225, tzxz =226, tzxt =227,
                     tzyx =228, tzyy =229, tzyz =230, tzyt =231,
                     tzzx =232, tzzy =233, tzzz =234, tzzt =235,
                     tztx =236, tzty =237, tztz =238, tztt =239,
                     ttxx =240, ttxy =241, ttxz =242, ttxt =243,
                     ttyx =244, ttyy =245, ttyz =246, ttyt =247,
                     ttzx =248, ttzy =249, ttzz =250, ttzt =251,
                     tttx =252, ttty =253, tttz =254, tttt =255 };
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
  template <class Ind>
  static constexpr Ind pascal_number(Ind&& i) { return 1; }
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
    if (list.size() > 1) {
      assert(list.size() == size() && "dimension size mismatch");
      size_t i = 0;
      for(auto it = list.begin(); it!=list.end(); ++it, ++i)
        data_[i]= *it;
    }
    else {
      auto it = list.begin();
      for(size_t i = 0; i < size(); ++i)
        data_[i]= *it;
      
    }
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

  // tensors of rank 3: access using Q[yxx], Q[zyx] etc. notation
  T& operator[](const enum tensor_indices::rank3_index ind) {
    assert (RANK == 3);
    using namespace tensor_indices;

    // generic case: no symmetry
    if constexpr (ST == symmetry_type::generic) {
      if constexpr (D == 1) {
        assert (ind == 0);
        return data_[0];
      }
      if constexpr (D == 2) {
        constexpr int remap[] = { 0, 4, -1, -1,
                                  2, 6, -1, -1,
                                  -1,-1,-1,-1,
                                  -1,-1,-1,-1,
                                  
                                  1, 5, -1, -1,
                                  3, 7};
        assert (ind < 22);
        assert (remap[ind] >= 0);
        return data_[remap[ind]];
      }
      if constexpr (D == 3) {
        constexpr int remap[] = { 0, 9,18, -1,
                                  3,12,21, -1,
                                  6,15,24,-1,
                                  -1,-1,-1,-1,
                                  
                                  1,10,19, -1,
                                  4,13,22, -1,
                                  7,16,25,-1,
                                  -1,-1,-1,-1,
                                  
                                  2,11,20, -1,
                                  5,14,23, -1,
                                  8,17,26};
        assert (ind < 43);
        assert (remap[ind] >= 0);
        return data_[remap[ind]];
      }
      if constexpr (D == 4) {
        constexpr int remap[] = { 0,16,32,48,
                                  4,20,36,52,
                                  8,24,40,56,
                                 12,28,44,60,
                                  
                                  1,17,33,49,
                                  5,21,37,53,
                                  9,25,41,57,
                                 13,29,45,61,
                                  
                                  2,18,34,50,
                                  6,22,38,54,
                                 10,26,42,58,
                                 14,30,46,62,
                                  
                                  3,19,35,51,
                                  7,23,39,55,
                                 11,27,43,59,
                                 15,31,47,63};
        assert (ind < 64);
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
                                  1, 2, -1, -1,
                                  -1,-1,-1,-1,
                                  -1,-1,-1,-1,
                                  
                                  1, 2, -1, -1,
                                  2, 3};
        assert (ind < 22);
        assert (remap[ind] >= 0);
        return data_[remap[ind]];
      }
      if constexpr (D == 3) {
        constexpr int remap[] = { 0, 1, 4, -1, 
                                  1, 2, 5, -1,
                                  4, 5, 7, -1,
                                  -1,-1,-1,-1,

                                  1, 2, 5, -1,
                                  2, 3, 6, -1,
                                  5, 6, 8, -1,
                                  -1,-1,-1,-1,

                                  4, 5, 7, -1,
                                  5, 6, 8, -1,
                                  7, 8, 9};
        assert (ind < 43);
        assert (remap[ind] >= 0);
        return data_[remap[ind]];
      }
      if constexpr (D == 4) {
        constexpr int remap[] = { 0, 1, 4, 10, 
                                  1, 2, 5, 11, 
                                  4, 5, 7, 13, 
                                  10, 11, 13, 16, 

                                  1, 2, 5, 11, 
                                  2, 3, 6, 12, 
                                  5, 6, 8, 14, 
                                  11, 12, 14, 17, 

                                  4, 5, 7, 13, 
                                  5, 6, 8, 14, 
                                  7, 8, 9, 15, 
                                  13, 14, 15, 18, 

                                  10, 11, 13, 16, 
                                  11, 12, 14, 17, 
                                  13, 14, 15, 18, 
                                  16, 17, 18, 19};

        assert (ind < 64);
        return data_[remap[ind]];
      }
    }

    assert (false);
  }

  // Pascal number: size()-th row, i-th diagonal
  template <class Ind>
  static constexpr Ind pascal_number(Ind&& i) {
    return i * tensor_u<T, ST, Ds...>::pascal_number(i + 1) / RANK;
  }

  // multi-index function: converts multi-index into flat index
  // constexpr for compile-time eval
  template <class Ind, class... Inds>
  static constexpr auto multiindex(Ind&& i, Inds&&... inds) {
    if constexpr (ST == symmetry_type::generic)
      return D * tensor_u<T, ST, Ds...>::multiindex(inds...) + i;

    if constexpr (ST == symmetry_type::symmetric) {
      if constexpr (RANK == 2)
        return symindx2(i, inds...);

      if constexpr (RANK == 3) 
        return symindx3(i, inds...);

      if constexpr (RANK == 4) 
        return symindx4(i, inds...);
    }

  }

  template <class Ind>
  static constexpr auto symindx2(Ind&& i, Ind&& j) {
    return (i>j) ? (i*(i+1)/2 + j) : (j*(j+1)/2 + i);
  }

  template <class Ind>
  static constexpr auto symindx3(Ind&& i, Ind&& j, Ind&& k) {
    return ((i>=j) && (i>=k)) ? (i*(i+1)*(i+2)/6 + symindx2(j,k)) :
           ((j>=i) && (j>=k)) ? (j*(j+1)*(j+2)/6 + symindx2(i,k)) :
                                (k*(k+1)*(k+2)/6 + symindx2(i,j));
  }

  template <class Ind>
  static constexpr auto symindx4(Ind&& i, Ind&& j, Ind&& k, Ind&& l) {
    return 
     ((i>=j)&&(i>=k)&&(i>=l)) ? (i*(i+1)*(i+2)*(i+3)/24 + symindx3(j,k,l)) :
     ((j>=i)&&(j>=k)&&(j>=l)) ? (j*(j+1)*(j+2)*(j+3)/24 + symindx3(i,k,l)) :
     ((k>=i)&&(k>=j)&&(k>=l)) ? (k*(k+1)*(k+2)*(k+3)/24 + symindx3(i,j,l)) :
                                (l*(l+1)*(l+2)*(l+3)/24 + symindx3(i,j,k));
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


/*!
  \function      operator == (tensor_u, tensor_u)
  \brief         Compare two tensors component-by-component

  \tparam T      Data type
  \tparam ST     Symmetry class
  \tparam Ds...  Dimensions of the vector space on which tensor is defined

  \param a       Tensor A
  \param b       Tensor B
 */
template <class T, symmetry_type ST, auto... Ds>
bool 
operator==(const tensor_u<T, ST, Ds...> & a, 
   const tensor_u<T, ST, Ds...> & b) {
  for(size_t i=0; i<tensor_u<T, ST, Ds...>::size(); ++i)
    if (a[i] != b[i]) return false;
  return true;
} // operator (==)


/*!
  \function      operator != (tensor_u, tensor_u)
  \brief         Compare two tensors with prejudice

  \tparam T      Data type
  \tparam ST     Symmetry class
  \tparam Ds...  Dimensions of the vector space on which tensor is defined

  \param a       Tensor A
  \param b       Tensor B
 */
template <class T, symmetry_type ST, auto... Ds>
bool 
operator!=(const tensor_u<T, ST, Ds...> & a, 
   const tensor_u<T, ST, Ds...> & b) {
  bool answer = false;
  for(size_t i=0; i<tensor_u<T, ST, Ds...>::size(); ++i)
    if (a[i] != b[i]) return true;
  return false;
} // operator (==)
} // namespace flecsi

#endif // TENSOR_H


/*
// Usage example
#include <iostream>
#include <iomanip>

int main() {

  using namespace std;
  using namespace flecsi;
  using namespace tensor_indices;

  cout << "--- Generic tensor: ---" << endl;
  using gen_tensor_t = tensor_u<double, symmetry_type::generic, 3, 3>;
  gen_tensor_t G{0};
  G[xy] = 1.2; G[yz] = 2.3; G[xz] = 1.3;
  G(1,1) = 2.0; G[xx] = 0.1;
  cout << "size of G: " << G.size() << endl;
  cout << "tensor G: "  << endl << G << endl;

  cout << "--- Symmetric tensor: ---" << endl;
  using sym_tensor_t = tensor_u<double, symmetry_type::symmetric, 3, 3>;
  sym_tensor_t S{0};
  S[xy] = 1.2; S[yz] = 2.3; S[xz] = 1.3;
  S(1,1) = 2.0; S[xx] = 0.1;
  cout << "size of S: " << S.size() << endl;
  cout << "tensor S: "  << endl << S << endl;
  cout << "accessing via the [xy.. etc.] operator: "  << endl
       << "S[xx]  S[xy]  S[xz]     ||"
       << setw(3) << S[xx] << " "
       << setw(3) << S[xy] << " "
       << setw(3) << S[xz] << " "
       << "||" << endl
       << "S[yx]  S[yy]  S[yz]  =  ||"
       << setw(3) << S[yx] << " "
       << setw(3) << S[yy] << " "
       << setw(3) << S[yz] << " "
       << "||" << endl
       << "S[zx]  S[zy]  S[zz]     ||"
       << setw(3) << S[zx] << " "
       << setw(3) << S[zy] << " "
       << setw(3) << S[zz] << " "
       << "||" << endl;
  cout << "accessing via the '(i,j)' operator:" << endl;
  for (int i=0; i<3; ++i) {
    cout << "S("<<i<<",0) S("<<i<<",1) S("<<i<<",2)";
    if (i==1) cout << "  =  ||"; else cout << "     ||";
    for (int j=0; j<3; ++j) cout << setw(3) << S(i,j) << " ";
    cout <<"||"<< endl;
  }

  cout << "--- Generic tensor, rank 3: ---" << endl;
  using gen2_tensor3_t = tensor_u<double, symmetry_type::generic, 2,2,2>;
  gen2_tensor3_t G3{0};
  G3[xxx]   = 111; G3[xxy]   = 112; G3[xyx]   = 121; G3[xyy]   = 122;
  G3(1,0,0) = 211; G3(1,0,1) = 212; G3(1,1,0) = 221; G3(1,1,1) = 222;
  cout << "size of G3: " << G3.size() << endl;
  cout << "tensor G3: "  << endl << G3 << endl;
  cout << "accessing via the [xyz.. etc.] operator: "  << endl;
  cout << "G3[xxx] = " << G3[xxx] << "; G3[xxy] = " << G3[xxy] << endl;
  cout << "G3[xyx] = " << G3[xyx] << "; G3[xyy] = " << G3[xyy] << endl;
  cout << "G3[yxx] = " << G3[yxx] << "; G3[yxy] = " << G3[yxy] << endl;
  cout << "G3[yyx] = " << G3[yyx] << "; G3[yyy] = " << G3[yyy] << endl;

  cout << "--- Generic tensor, rank 3, dimension 3: ---" << endl;
  using gen3_tensor3_t = tensor_u<double, symmetry_type::generic, 3,3,3>;
  gen3_tensor3_t Q3{0};
  Q3[xxx]   = 111; Q3[xxy]   = 112; Q3[xyx]   = 121; Q3[xyy]   = 122;
  Q3(1,0,0) = 211; Q3(1,0,1) = 212; Q3(1,1,0) = 221; Q3(1,1,1) = 222;
  cout << "size of Q3: " << Q3.size() << endl;
  cout << "tensor Q3: "  << endl << Q3 << endl;
  cout << "accessing via the [xyz.. etc.] operator: "  << endl;
  cout << "Q3[xxx] = " << Q3[xxx] << "; Q3[xxy] = " << Q3[xxy] << endl;
  cout << "Q3[xyx] = " << Q3[xyx] << "; Q3[xyy] = " << Q3[xyy] << endl;
  cout << "Q3[yxx] = " << Q3[yxx] << "; Q3[yxy] = " << Q3[yxy] << endl;
  cout << "Q3[yyx] = " << Q3[yyx] << "; Q3[yyy] = " << Q3[yyy] << endl;

  cout << "--- Symmetric tensor, rank 3: ---" << endl;
  using sym2_tensor3_t = tensor_u<double, symmetry_type::symmetric, 2,2,2>;
  sym2_tensor3_t S3{0};
  S3[xxx]   = 111; S3[xxy]   = 112; S3[xyy]   = 122; 
  S3(1,1,1) = 222;
  cout << "size of S3: " << S3.size() << endl;
  cout << "tensor S3: "  << endl << S3 << endl;
  cout << "accessing via the [xyz.. etc.] operator: "  << endl;
  cout << "S3[xxx] = " << S3[xxx] << "; S3[xxy] = " << S3[xxy] << endl;
  cout << "S3[xyx] = " << S3[xyx] << "; S3[xyy] = " << S3[xyy] << endl;
  cout << "S3[yxx] = " << S3[yxx] << "; S3[yxy] = " << S3[yxy] << endl;
  cout << "S3[yyx] = " << S3[yyx] << "; S3[yyy] = " << S3[yyy] << endl;

  cout << "--- Symmetric tensor, rank 3, dimension 3: ---" << endl;
  using sym3_tensor3_t = tensor_u<double, symmetry_type::symmetric,3,3,3>;
  sym3_tensor3_t Z3{0};
  Z3[xxx]   = 111; Z3(_y,_y,_y) = 222; Z3(2,2,2) = 333;
  Z3[xxy]   = 112; Z3[xyy]   = 122; Z3[xxz]   = 113; Z3[xzz]  = 133;
  Z3[xyz]   = 123; Z3[yyz]   = 223; Z3[yzz]   = 233; 
  cout << "size of Z3: " << Z3.size() << endl;
  cout << "tensor Z3: "  << endl << Z3 << endl;
  cout << "accessing via the [xyz.. etc.] operator: "  << endl;
  cout << "Z3[xxx] = " << Z3[xxx] 
       << "; Z3[xxy] = " << Z3[xxy]
       << "; Z3[xxz] = " << Z3[xxz] << endl;
  cout << "Z3[xyx] = " << Z3[xyx] 
       << "; Z3[xyy] = " << Z3[xyy]
       << "; Z3[xyz] = " << Z3[xyz] << endl;
  cout << "Z3[xzx] = " << Z3[xzx] 
       << "; Z3[xzy] = " << Z3[xzy]
       << "; Z3[xzz] = " << Z3[xzz] << endl;
  cout << endl;
  cout << "Z3[yxx] = " << Z3[yxx] 
       << "; Z3[yxy] = " << Z3[yxy]
       << "; Z3[yxz] = " << Z3[yxz] << endl;
  cout << "Z3[yyx] = " << Z3[yyx] 
       << "; Z3[yyy] = " << Z3[yyy]
       << "; Z3[yyz] = " << Z3[yyz] << endl;
  cout << "Z3[yzx] = " << Z3[yzx] 
       << "; Z3[yzy] = " << Z3[yzy]
       << "; Z3[yzz] = " << Z3[yzz] << endl;
  cout << endl;
  cout << "Z3[zxx] = " << Z3[zxx] 
       << "; Z3[zxy] = " << Z3[zxy]
       << "; Z3[zxz] = " << Z3[zxz] << endl;
  cout << "Z3[zyx] = " << Z3[zyx] 
       << "; Z3[zyy] = " << Z3[zyy]
       << "; Z3[zyz] = " << Z3[zyz] << endl;
  cout << "Z3[zzx] = " << Z3[zzx] 
       << "; Z3[zzy] = " << Z3[zzy]
       << "; Z3[zzz] = " << Z3[zzz] << endl;

  cout << "--- Symmetric tensor, rank 3, dimension 4: ---" << endl;
  using sym3_tensor4_t = tensor_u<double, symmetry_type::symmetric,4,4,4>;
  for (int i=0;i<4; ++i) {
    for (int j=0;j<4; ++j) {
      for (int k=0;k<4; ++k) {
         cout << sym3_tensor4_t::multiindex(i,j,k) << ", ";
      }
      cout << endl;
    }
    cout << endl;
  }
  
  //cout << "--- Higher-rank symmetric tensor: ---" << endl;
  //using sym_tensor4_t = tensor_u<double, symmetry_type::symmetric, 3,3,3,3>;
  //sym_tensor4_t Q3{0};
  //cout << "--- Q3 size = " << Q3.size() << endl;
  //Q3(0,0,0,0) = 3333;
  //Q3(0,1,0,1) = 12;
  //Q3(1,2,1,0) = 31;
  //for (int m=0;m<3;++m) {
  //  for (int i=0;i<3; ++i) {
  //    for (int j=0;j<3; ++j) {
  //      for (int k=0;k<3; ++k) {
  //         cout << Q3(m,i,j,k) << " ";
  //      }
  //      cout << endl;
  //    }
  //    cout << " ----------- " << endl;
  //  }
  //  cout << "=============" << endl;
  //}
}
*/

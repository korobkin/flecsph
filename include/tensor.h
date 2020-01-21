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

    if constexpr (ST == symmetry_type::symmetric)
      return pascal_number(D);

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


  // tensors of rank 4: access using Q[ytxx], Q[zytx] etc. notation
  T& operator[](const enum tensor_indices::rank4_index ind) {
    assert (RANK == 4);
    using namespace tensor_indices;

    // generic case: no symmetry
    if constexpr (ST == symmetry_type::generic) {
      if constexpr (D == 1) {
        assert (ind == 0);
        return data_[0];
      }
      if constexpr (D == 2) {
        constexpr int remap[] = {  0,  8, -1, -1,
                                   4, 12, -1, -1,
                                  -1, -1, -1, -1,
                                  -1, -1, -1, -1,

                                   2, 10, -1, -1,
                                   6, 14, -1, -1,
                                  -1, -1, -1, -1,
                                  -1, -1, -1, -1,

                                  -1, -1, -1, -1,
                                  -1, -1, -1, -1,
                                  -1, -1, -1, -1,
                                  -1, -1, -1, -1,

                                  -1, -1, -1, -1,
                                  -1, -1, -1, -1,
                                  -1, -1, -1, -1,
                                  -1, -1, -1, -1,

                                   1,  9, -1, -1,
                                   5, 13, -1, -1,
                                  -1, -1, -1, -1,
                                  -1, -1, -1, -1,

                                  3, 11, -1, -1,
                                  7, 15};        
        assert (ind < 86);
        assert (remap[ind] >= 0);
        return data_[remap[ind]];
      }
      if constexpr (D == 3) {
        constexpr int remap[] = {  0, 27, 54, -1,
                                   9, 36, 63, -1,
                                  18, 45, 72, -1,
                                  -1, -1, -1, -1,

                                   3, 30, 57, -1,
                                  12, 39, 66, -1,
                                  21, 48, 75, -1,
                                  -1, -1, -1, -1,

                                   6, 33, 60, -1,
                                  15, 42, 69, -1,
                                  24, 51, 78, -1,
                                  -1, -1, -1, -1,

                                  -1, -1, -1, -1,
                                  -1, -1, -1, -1,
                                  -1, -1, -1, -1,
                                  -1, -1, -1, -1,

                                   1, 28, 55, -1,
                                  10, 37, 64, -1,
                                  19, 46, 73, -1,
                                  -1, -1, -1, -1,

                                   4, 31, 58, -1,
                                  13, 40, 67, -1,
                                  22, 49, 76, -1,
                                  -1, -1, -1, -1,

                                   7, 34, 61, -1,
                                  16, 43, 70, -1,
                                  25, 52, 79, -1,
                                  -1, -1, -1, -1,

                                  -1, -1, -1, -1,
                                  -1, -1, -1, -1,
                                  -1, -1, -1, -1,
                                  -1, -1, -1, -1,

                                   2, 29, 56, -1,
                                  11, 38, 65, -1,
                                  20, 47, 74, -1,
                                  -1, -1, -1, -1,

                                   5, 32, 59, -1,
                                  14, 41, 68, -1,
                                  23, 50, 77, -1,
                                  -1, -1, -1, -1,

                                   8, 35, 62, -1,
                                  17, 44, 71, -1,
                                  26, 53, 80};
        assert (ind < 171);
        assert (remap[ind] >= 0);
        return data_[remap[ind]];
      }
      if constexpr (D == 4) {
        constexpr int remap[] = { 0, 64, 128, 192, 
                                 16, 80, 144, 208, 
                                 32, 96, 160, 224, 
                                 48,112, 176, 240, 

                                  4, 68, 132, 196, 
                                 20, 84, 148, 212, 
                                 36,100, 164, 228, 
                                 52,116, 180, 244, 

                                  8, 72, 136, 200, 
                                 24, 88, 152, 216, 
                                 40,104, 168, 232, 
                                 56,120, 184, 248, 

                                 12, 76, 140, 204, 
                                 28, 92, 156, 220, 
                                 44,108, 172, 236, 
                                 60,124, 188, 252, 


                                  1, 65, 129, 193, 
                                 17, 81, 145, 209, 
                                 33, 97, 161, 225, 
                                 49,113, 177, 241, 

                                  5, 69, 133, 197, 
                                 21, 85, 149, 213, 
                                 37,101, 165, 229, 
                                 53,117, 181, 245, 

                                  9, 73, 137, 201, 
                                 25, 89, 153, 217, 
                                 41,105, 169, 233, 
                                 57,121, 185, 249, 

                                 13, 77, 141, 205, 
                                 29, 93, 157, 221, 
                                 45,109, 173, 237, 
                                 61,125, 189, 253, 


                                  2, 66, 130, 194, 
                                 18, 82, 146, 210, 
                                 34, 98, 162, 226, 
                                 50,114, 178, 242, 

                                  6, 70, 134, 198, 
                                 22, 86, 150, 214, 
                                 38,102, 166, 230, 
                                 54,118, 182, 246, 

                                 10, 74, 138, 202, 
                                 26, 90, 154, 218, 
                                 42,106, 170, 234, 
                                 58,122, 186, 250, 

                                 14, 78, 142, 206, 
                                 30, 94, 158, 222, 
                                 46,110, 174, 238, 
                                 62,126, 190, 254, 


                                  3, 67, 131, 195, 
                                 19, 83, 147, 211, 
                                 35, 99, 163, 227, 
                                 51,115, 179, 243, 

                                  7, 71, 135, 199, 
                                 23, 87, 151, 215, 
                                 39,103, 167, 231, 
                                 55,119, 183, 247, 

                                 11, 75, 139, 203, 
                                 27, 91, 155, 219, 
                                 43,107, 171, 235, 
                                 59,123, 187, 251, 

                                 15, 79, 143, 207, 
                                 31, 95, 159, 223, 
                                 47,111, 175, 239, 
                                 63,127, 191, 255};
        assert (ind < 256);
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
        constexpr int remap[] = {   0, 1, -1, -1,
                                    1, 2, -1, -1,
                                  -1, -1, -1, -1,
                                  -1, -1, -1, -1,

                                   1,  2, -1, -1,
                                   2,  3, -1, -1,
                                  -1, -1, -1, -1,
                                  -1, -1, -1, -1,

                                  -1, -1, -1, -1,
                                  -1, -1, -1, -1,
                                  -1, -1, -1, -1,
                                  -1, -1, -1, -1,

                                  -1, -1, -1, -1,
                                  -1, -1, -1, -1,
                                  -1, -1, -1, -1,
                                  -1, -1, -1, -1,

                                   1,  2, -1, -1,
                                   2,  3, -1, -1,
                                  -1, -1, -1, -1,
                                  -1, -1, -1, -1,

                                   2,  3, -1, -1,
                                   3, 4};
        assert (ind < 86);
        assert (remap[ind] >= 0);
        return data_[remap[ind]];
      }
      if constexpr (D == 3) {
        constexpr int remap[] = {   0,  1,  5, -1,
                                    1,  2,  6, -1,
                                    5,  6,  9, -1,
                                   -1, -1, -1, -1,

                                    1,  2,  6, -1,
                                    2,  3,  7, -1,
                                    6,  7, 10, -1,
                                   -1, -1, -1, -1,

                                    5,  6,  9, -1,
                                    6,  7, 10, -1,
                                    9, 10, 12, -1,
                                   -1, -1, -1, -1,

                                   -1, -1, -1, -1,
                                   -1, -1, -1, -1,
                                   -1, -1, -1, -1,
                                   -1, -1, -1, -1,

                                    1,  2,  6, -1,
                                    2,  3,  7, -1,
                                    6,  7, 10, -1,
                                   -1, -1, -1, -1,

                                    2,  3,  7, -1,
                                    3,  4,  8, -1,
                                    7,  8, 11, -1,
                                   -1, -1, -1, -1,

                                    6,  7, 10, -1,
                                    7,  8, 11, -1,
                                   10, 11, 13, -1,
                                   -1, -1, -1, -1,

                                   -1, -1, -1, -1,
                                   -1, -1, -1, -1,
                                   -1, -1, -1, -1,
                                   -1, -1, -1, -1,

                                    5,  6,  9, -1,
                                    6,  7, 10, -1,
                                    9, 10, 12, -1,
                                   -1, -1, -1, -1,

                                    6,  7, 10, -1,
                                    7,  8, 11, -1,
                                   10, 11, 13, -1,
                                   -1, -1, -1, -1,

                                    9, 10, 12, -1,
                                   10, 11, 13, -1,
                                   12, 13, 14};
        assert (ind < 171);
        assert (remap[ind] >= 0);
        return data_[remap[ind]];
      }
      if constexpr (D == 4) {
        constexpr int remap[] = {   0,  1,  5, 15, 
                                    1,  2,  6, 16, 
                                    5,  6,  9, 19, 
                                   15, 16, 19, 25, 

                                    1,  2,  6, 16, 
                                    2,  3,  7, 17, 
                                    6,  7, 10, 20, 
                                   16, 17, 20, 26, 

                                    5,  6,  9, 19, 
                                    6,  7, 10, 20, 
                                    9, 10, 12, 22, 
                                   19, 20, 22, 28, 

                                   15, 16, 19, 25, 
                                   16, 17, 20, 26, 
                                   19, 20, 22, 28, 
                                   25, 26, 28, 31, 


                                    1,  2,  6, 16, 
                                    2,  3,  7, 17, 
                                    6,  7, 10, 20, 
                                   16, 17, 20, 26, 

                                    2,  3,  7, 17, 
                                    3,  4,  8, 18, 
                                    7,  8, 11, 21, 
                                   17, 18, 21, 27, 

                                    6,  7, 10, 20, 
                                    7,  8, 11, 21, 
                                   10, 11, 13, 23, 
                                   20, 21, 23, 29, 

                                   16, 17, 20, 26, 
                                   17, 18, 21, 27, 
                                   20, 21, 23, 29, 
                                   26, 27, 29, 32, 


                                    5,  6,  9, 19, 
                                    6,  7, 10, 20, 
                                    9, 10, 12, 22, 
                                   19, 20, 22, 28, 

                                    6,  7, 10, 20, 
                                    7,  8, 11, 21, 
                                   10, 11, 13, 23, 
                                   20, 21, 23, 29, 

                                    9, 10, 12, 22, 
                                   10, 11, 13, 23, 
                                   12, 13, 14, 24, 
                                   22, 23, 24, 30, 

                                   19, 20, 22, 28, 
                                   20, 21, 23, 29, 
                                   22, 23, 24, 30, 
                                   28, 29, 30, 33, 


                                   15, 16, 19, 25, 
                                   16, 17, 20, 26, 
                                   19, 20, 22, 28, 
                                   25, 26, 28, 31, 

                                   16, 17, 20, 26, 
                                   17, 18, 21, 27, 
                                   20, 21, 23, 29, 
                                   26, 27, 29, 32, 

                                   19, 20, 22, 28, 
                                   20, 21, 23, 29, 
                                   22, 23, 24, 30, 
                                   28, 29, 30, 33, 

                                   25, 26, 28, 31, 
                                   26, 27, 29, 32, 
                                   28, 29, 30, 33, 
                                   31, 32, 33, 34};

        assert (ind < 256);
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

#define def_operator(op)                                                       \
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

#define def_operator_type(op)                                                  \
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

  def_operator(+=);

  //--------------------------------------------------------------------------//
  //! Addition/Assignment operator.
  //--------------------------------------------------------------------------//

  def_operator_type(+=);

  //--------------------------------------------------------------------------//
  //! Subtraction/Assignment operator.
  //--------------------------------------------------------------------------//

  def_operator(-=);

  //--------------------------------------------------------------------------//
  //! Subtraction/Assignment operator.
  //--------------------------------------------------------------------------//

  def_operator_type(-=);

  //--------------------------------------------------------------------------//
  //! Division/Assignment operator.
  //--------------------------------------------------------------------------//

  def_operator_type(*=);

  //--------------------------------------------------------------------------//
  //! Division/Assignment operator.
  //--------------------------------------------------------------------------//

  def_operator_type(/=);

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

  cout << "--- Generic tensor, rank 4, dimension 2: ---" << endl;
  using gen2_tensor4_t = tensor_u<double, symmetry_type::generic,2,2,2,2>;
  gen2_tensor4_t Q4{0};
  Q4[xxxx]=1111; Q4[xxxy]=1112; Q4[xxyx]=1121; Q4[xxyy]=1122;
  Q4(0,1,0,0)=1211; Q4(0,1,0,1)=1212; Q4(0,1,1,0)=1221; Q4(0,1,1,1)=1222;
  Q4(1,0,0,0)=2111; Q4(1,0,0,1)=2112; Q4(1,0,1,0)=2121; Q4(1,0,1,1)=2122;
  Q4(1,1,0,0)=2211; Q4(1,1,0,1)=2212; Q4(1,1,1,0)=2221; Q4(1,1,1,1)=2222;
  cout << "size of Q4: " << Q4.size() << endl;
  cout << "tensor Q4: "  << endl << Q4 << endl;
  cout << "accessing via the [xyxx.. etc.] operator: "  << endl;
  cout << "Q4[xxxx] = " << Q4[xxxx] << "; Q4[xxxy] = " << Q4[xxxy] << "; "
       << "Q4[xxyx] = " << Q4[xxyx] << "; Q4[xxyy] = " << Q4[xxyy] << endl;
  cout << "Q4[xyxx] = " << Q4[xyxx] << "; Q4[xyxy] = " << Q4[xyxy] << "; "
       << "Q4[xyyx] = " << Q4[xyyx] << "; Q4[xyyy] = " << Q4[xyyy] << endl;
  cout << "Q4[yxxx] = " << Q4[yxxx] << "; Q4[yxxy] = " << Q4[yxxy] << "; "
       << "Q4[yxyx] = " << Q4[yxyx] << "; Q4[yxyy] = " << Q4[yxyy] << endl;
  cout << "Q4[yyxx] = " << Q4[yyxx] << "; Q4[yyxy] = " << Q4[yyxy] << "; "
       << "Q4[yyyx] = " << Q4[yyyx] << "; Q4[yyyy] = " << Q4[yyyy] << endl;

  cout << "--- Generic tensor, rank 4, dimension 3: ---" << endl;
  using gen3_tensor4_t = tensor_u<double, symmetry_type::generic,3,3,3,3>;
  gen3_tensor4_t R4{0};
  R4[xxxx]=1111; R4[xxxy]=1112; R4[xxxz]=1113; 
  R4[xxyx]=1121; R4[xxyy]=1122; R4[xxyz]=1123; 
  R4[xxzx]=1131; R4[xxzy]=1132; R4[xxzz]=1133; 

  R4[xyxx]=1211; R4[xyxy]=1212; R4[xyxz]=1213; 
  R4[xyyx]=1221; R4[xyyy]=1222; R4[xyyz]=1223; 
  R4[xyzx]=1231; R4[xyzy]=1232; R4[xyzz]=1233; 

  R4[xzxx]=1311; R4[xzxy]=1312; R4[xzxz]=1313; 
  R4[xzyx]=1321; R4[xzyy]=1322; R4[xzyz]=1323; 
  R4[xzzx]=1331; R4[xzzy]=1332; R4[xzzz]=1333; 


  R4(1,0,0,0)=2111; R4(1,0,0,1)=2112; R4(1,0,0,2)=2113; 
  R4(1,0,1,0)=2121; R4(1,0,1,1)=2122; R4(1,0,1,2)=2123; 
  R4(1,0,2,0)=2131; R4(1,0,2,1)=2132; R4(1,0,2,2)=2133; 

  R4[yyxx]=2211; R4[yyxy]=2212; R4[yyxz]=2213; 
  R4[yyyx]=2221; R4[yyyy]=2222; R4[yyyz]=2223; 
  R4[yyzx]=2231; R4[yyzy]=2232; R4[yyzz]=2233; 

  R4[yzxx]=2311; R4[yzxy]=2312; R4[yzxz]=2313; 
  R4[yzyx]=2321; R4[yzyy]=2322; R4[yzyz]=2323; 
  R4[yzzx]=2331; R4[yzzy]=2332; R4[yzzz]=2333; 


  R4[zxxx]=3111; R4[zxxy]=3112; R4[zxxz]=3113; 
  R4[zxyx]=3121; R4[zxyy]=3122; R4[zxyz]=3123; 
  R4[zxzx]=3131; R4[zxzy]=3132; R4[zxzz]=3133; 

  R4[zyxx]=3211; R4[zyxy]=3212; R4[zyxz]=3213; 
  R4[zyyx]=3221; R4[zyyy]=3222; R4[zyyz]=3223; 
  R4[zyzx]=3231; R4[zyzy]=3232; R4[zyzz]=3233; 

  R4[zzxx]=3311; R4[zzxy]=3312; R4[zzxz]=3313; 
  R4[zzyx]=3321; R4[zzyy]=3322; R4[zzyz]=3323; 
  R4[zzzx]=3331; R4[zzzy]=3332; R4[zzzz]=3333; 

  cout << "size of R4: " << R4.size() << endl;
  cout << "tensor R4: "  << endl << R4 << endl;
  cout << "accessing via the [xyxx.. etc.] operator: "  << endl;
  cout << "R4[xxxx] = " << R4[xxxx] << "; "
       << "R4[xxxy] = " << R4[xxxy] << "; "
       << "R4[xxxz] = " << R4[xxxz] << endl
       << "R4[xxyx] = " << R4[xxyx] << "; "
       << "R4[xxyy] = " << R4[xxyy] << "; "
       << "R4[xxyz] = " << R4[xxyz] << endl
       << "R4[xxzx] = " << R4[xxzx] << "; "
       << "R4[xxzy] = " << R4[xxzy] << "; "
       << "R4[xxzz] = " << R4[xxzz] << endl
       << endl;
  cout << "R4[xyxx] = " << R4[xyxx] << "; "
       << "R4[xyxy] = " << R4[xyxy] << "; "
       << "R4[xyxz] = " << R4[xyxz] << endl
       << "R4[xyyx] = " << R4[xyyx] << "; "
       << "R4[xyyy] = " << R4[xyyy] << "; "
       << "R4[xyyz] = " << R4[xyyz] << endl
       << "R4[xyzx] = " << R4[xyzx] << "; "
       << "R4[xyzy] = " << R4[xyzy] << "; "
       << "R4[xyzz] = " << R4[xyzz] << endl
       << endl;
  cout << "R4[xzxx] = " << R4[xzxx] << "; "
       << "R4[xzxy] = " << R4[xzxy] << "; "
       << "R4[xzxz] = " << R4[xzxz] << endl
       << "R4[xzyx] = " << R4[xzyx] << "; "
       << "R4[xzyy] = " << R4[xzyy] << "; "
       << "R4[xzyz] = " << R4[xzyz] << endl
       << "R4[xzzx] = " << R4[xzzx] << "; "
       << "R4[xzzy] = " << R4[xzzy] << "; "
       << "R4[xzzz] = " << R4[xzzz] << endl
       << endl;
  cout << endl;

  cout << "R4[yxxx] = " << R4[yxxx] << "; "
       << "R4[yxxy] = " << R4[yxxy] << "; "
       << "R4[yxxz] = " << R4[yxxz] << endl
       << "R4[yxyx] = " << R4[yxyx] << "; "
       << "R4[yxyy] = " << R4[yxyy] << "; "
       << "R4[yxyz] = " << R4[yxyz] << endl
       << "R4[yxzx] = " << R4[yxzx] << "; "
       << "R4[yxzy] = " << R4[yxzy] << "; "
       << "R4[yxzz] = " << R4[yxzz] << endl
       << endl;
  cout << "R4[yyxx] = " << R4[yyxx] << "; "
       << "R4[yyxy] = " << R4[yyxy] << "; "
       << "R4[yyxz] = " << R4[yyxz] << endl
       << "R4[yyyx] = " << R4[yyyx] << "; "
       << "R4[yyyy] = " << R4[yyyy] << "; "
       << "R4[yyyz] = " << R4[yyyz] << endl
       << "R4[yyzx] = " << R4[yyzx] << "; "
       << "R4[yyzy] = " << R4[yyzy] << "; "
       << "R4[yyzz] = " << R4[yyzz] << endl
       << endl;
  cout << "R4[yzxx] = " << R4[yzxx] << "; "
       << "R4[yzxy] = " << R4[yzxy] << "; "
       << "R4[yzxz] = " << R4[yzxz] << endl
       << "R4[yzyx] = " << R4[yzyx] << "; "
       << "R4[yzyy] = " << R4[yzyy] << "; "
       << "R4[yzyz] = " << R4[yzyz] << endl
       << "R4[yzzx] = " << R4[yzzx] << "; "
       << "R4[yzzy] = " << R4[yzzy] << "; "
       << "R4[yzzz] = " << R4[yzzz] << endl
       << endl;
  cout << endl;

  cout << "R4[zxxx] = " << R4[zxxx] << "; "
       << "R4[zxxy] = " << R4[zxxy] << "; "
       << "R4[zxxz] = " << R4[zxxz] << endl
       << "R4[zxyx] = " << R4[zxyx] << "; "
       << "R4[zxyy] = " << R4[zxyy] << "; "
       << "R4[zxyz] = " << R4[zxyz] << endl
       << "R4[zxzx] = " << R4[zxzx] << "; "
       << "R4[zxzy] = " << R4[zxzy] << "; "
       << "R4[zxzz] = " << R4[zxzz] << endl
       << endl;
  cout << "R4[zyxx] = " << R4[zyxx] << "; "
       << "R4[zyxy] = " << R4[zyxy] << "; "
       << "R4[zyxz] = " << R4[zyxz] << endl
       << "R4[zyyx] = " << R4[zyyx] << "; "
       << "R4[zyyy] = " << R4[zyyy] << "; "
       << "R4[zyyz] = " << R4[zyyz] << endl
       << "R4[zyzx] = " << R4[zyzx] << "; "
       << "R4[zyzy] = " << R4[zyzy] << "; "
       << "R4[zyzz] = " << R4[zyzz] << endl
       << endl;
  cout << "R4[zzxx] = " << R4[zzxx] << "; "
       << "R4[zzxy] = " << R4[zzxy] << "; "
       << "R4[zzxz] = " << R4[zzxz] << endl
       << "R4[zzyx] = " << R4[zzyx] << "; "
       << "R4[zzyy] = " << R4[zzyy] << "; "
       << "R4[zzyz] = " << R4[zzyz] << endl
       << "R4[zzzx] = " << R4[zzzx] << "; "
       << "R4[zzzy] = " << R4[zzzy] << "; "
       << "R4[zzzz] = " << R4[zzzz] << endl
       << endl;

  cout << "--- Symmetric tensor, rank 4, dimension 2: ---" << endl;
  using sym2_tensor4_t = tensor_u<double, symmetry_type::symmetric, 2,2,2,2>;
  sym2_tensor4_t S4{0};
  S4[xxxx] = 1111; S4[yyyy] = 2222;
  S4(0,0,0,1) = 1112; 
  S4(0,0,1,1) = 1122; 
  S4(1,0,1,1) = 1222; 
  cout << "size of S4: " << S4.size() << endl;
  cout << "tensor S4: "  << endl << S4 << endl;
  cout << "accessing via the [xyxx.. etc.] operator: "  << endl;
  cout << "S4[xxxx] = " << S4[xxxx] << "; S4[xxxy] = " << S4[xxxy] << "; "
       << "S4[xxyx] = " << S4[xxyx] << "; S4[xxyy] = " << S4[xxyy] << endl;
  cout << "S4[xyxx] = " << S4[xyxx] << "; S4[xyxy] = " << S4[xyxy] << "; "
       << "S4[xyyx] = " << S4[xyyx] << "; S4[xyyy] = " << S4[xyyy] << endl;
  cout << "S4[yxxx] = " << S4[yxxx] << "; S4[yxxy] = " << S4[yxxy] << "; "
       << "S4[yxyx] = " << S4[yxyx] << "; S4[yxyy] = " << S4[yxyy] << endl;
  cout << "S4[yyxx] = " << S4[yyxx] << "; S4[yyxy] = " << S4[yyxy] << "; "
       << "S4[yyyx] = " << S4[yyyx] << "; S4[yyyy] = " << S4[yyyy] << endl;

  cout << "--- Generic tensor, rank 4, dimension 3: ---" << endl;
  using sym3_tensor4_t = tensor_u<double, symmetry_type::symmetric,3,3,3,3>;
  sym3_tensor4_t Z4{0};
  Z4[xxxx]=1111, Z4[yyyy]=2222; Z4[zzzz]=3333;
  Z4[xxxy]=1112, Z4[xxyy]=1122; Z4[xyyy]=1222;
  Z4[xxxz]=1113, Z4[xxzz]=1133; Z4[xzzz]=1333;
  Z4[yyyz]=2223, Z4[yyzz]=2233; Z4[yzzz]=2333;
  Z4[xyzz]=1233, Z4[xyyz]=1223; Z4[xxyz]=1123;

  cout << "size of Z4: " << Z4.size() << endl;
  cout << "tensor Z4: "  << endl << Z4 << endl;
  cout << "accessing via the [xyxx.. etc.] operator: "  << endl;
  cout << "Z4[xxxx] = " << Z4[xxxx] << "; "
       << "Z4[xxxy] = " << Z4[xxxy] << "; "
       << "Z4[xxxz] = " << Z4[xxxz] << endl
       << "Z4[xxyx] = " << Z4[xxyx] << "; "
       << "Z4[xxyy] = " << Z4[xxyy] << "; "
       << "Z4[xxyz] = " << Z4[xxyz] << endl
       << "Z4[xxzx] = " << Z4[xxzx] << "; "
       << "Z4[xxzy] = " << Z4[xxzy] << "; "
       << "Z4[xxzz] = " << Z4[xxzz] << endl
       << endl;
  cout << "Z4[xyxx] = " << Z4[xyxx] << "; "
       << "Z4[xyxy] = " << Z4[xyxy] << "; "
       << "Z4[xyxz] = " << Z4[xyxz] << endl
       << "Z4[xyyx] = " << Z4[xyyx] << "; "
       << "Z4[xyyy] = " << Z4[xyyy] << "; "
       << "Z4[xyyz] = " << Z4[xyyz] << endl
       << "Z4[xyzx] = " << Z4[xyzx] << "; "
       << "Z4[xyzy] = " << Z4[xyzy] << "; "
       << "Z4[xyzz] = " << Z4[xyzz] << endl
       << endl;
  cout << "Z4[xzxx] = " << Z4[xzxx] << "; "
       << "Z4[xzxy] = " << Z4[xzxy] << "; "
       << "Z4[xzxz] = " << Z4[xzxz] << endl
       << "Z4[xzyx] = " << Z4[xzyx] << "; "
       << "Z4[xzyy] = " << Z4[xzyy] << "; "
       << "Z4[xzyz] = " << Z4[xzyz] << endl
       << "Z4[xzzx] = " << Z4[xzzx] << "; "
       << "Z4[xzzy] = " << Z4[xzzy] << "; "
       << "Z4[xzzz] = " << Z4[xzzz] << endl
       << endl;
  cout << endl;

  cout << "Z4[yxxx] = " << Z4[yxxx] << "; "
       << "Z4[yxxy] = " << Z4[yxxy] << "; "
       << "Z4[yxxz] = " << Z4[yxxz] << endl
       << "Z4[yxyx] = " << Z4[yxyx] << "; "
       << "Z4[yxyy] = " << Z4[yxyy] << "; "
       << "Z4[yxyz] = " << Z4[yxyz] << endl
       << "Z4[yxzx] = " << Z4[yxzx] << "; "
       << "Z4[yxzy] = " << Z4[yxzy] << "; "
       << "Z4[yxzz] = " << Z4[yxzz] << endl
       << endl;
  cout << "Z4[yyxx] = " << Z4[yyxx] << "; "
       << "Z4[yyxy] = " << Z4[yyxy] << "; "
       << "Z4[yyxz] = " << Z4[yyxz] << endl
       << "Z4[yyyx] = " << Z4[yyyx] << "; "
       << "Z4[yyyy] = " << Z4[yyyy] << "; "
       << "Z4[yyyz] = " << Z4[yyyz] << endl
       << "Z4[yyzx] = " << Z4[yyzx] << "; "
       << "Z4[yyzy] = " << Z4[yyzy] << "; "
       << "Z4[yyzz] = " << Z4[yyzz] << endl
       << endl;
  cout << "Z4[yzxx] = " << Z4[yzxx] << "; "
       << "Z4[yzxy] = " << Z4[yzxy] << "; "
       << "Z4[yzxz] = " << Z4[yzxz] << endl
       << "Z4[yzyx] = " << Z4[yzyx] << "; "
       << "Z4[yzyy] = " << Z4[yzyy] << "; "
       << "Z4[yzyz] = " << Z4[yzyz] << endl
       << "Z4[yzzx] = " << Z4[yzzx] << "; "
       << "Z4[yzzy] = " << Z4[yzzy] << "; "
       << "Z4[yzzz] = " << Z4[yzzz] << endl
       << endl;
  cout << endl;

  cout << "Z4[zxxx] = " << Z4[zxxx] << "; "
       << "Z4[zxxy] = " << Z4[zxxy] << "; "
       << "Z4[zxxz] = " << Z4[zxxz] << endl
       << "Z4[zxyx] = " << Z4[zxyx] << "; "
       << "Z4[zxyy] = " << Z4[zxyy] << "; "
       << "Z4[zxyz] = " << Z4[zxyz] << endl
       << "Z4[zxzx] = " << Z4[zxzx] << "; "
       << "Z4[zxzy] = " << Z4[zxzy] << "; "
       << "Z4[zxzz] = " << Z4[zxzz] << endl
       << endl;
  cout << "Z4[zyxx] = " << Z4[zyxx] << "; "
       << "Z4[zyxy] = " << Z4[zyxy] << "; "
       << "Z4[zyxz] = " << Z4[zyxz] << endl
       << "Z4[zyyx] = " << Z4[zyyx] << "; "
       << "Z4[zyyy] = " << Z4[zyyy] << "; "
       << "Z4[zyyz] = " << Z4[zyyz] << endl
       << "Z4[zyzx] = " << Z4[zyzx] << "; "
       << "Z4[zyzy] = " << Z4[zyzy] << "; "
       << "Z4[zyzz] = " << Z4[zyzz] << endl
       << endl;
  cout << "Z4[zzxx] = " << Z4[zzxx] << "; "
       << "Z4[zzxy] = " << Z4[zzxy] << "; "
       << "Z4[zzxz] = " << Z4[zzxz] << endl
       << "Z4[zzyx] = " << Z4[zzyx] << "; "
       << "Z4[zzyy] = " << Z4[zzyy] << "; "
       << "Z4[zzyz] = " << Z4[zzyz] << endl
       << "Z4[zzzx] = " << Z4[zzzx] << "; "
       << "Z4[zzzy] = " << Z4[zzzy] << "; "
       << "Z4[zzzz] = " << Z4[zzzz] << endl
       << endl;

  //cout << "--- Symmetric tensor, rank 3, dimension 4: ---" << endl;
  //using sym4_tensor3_t = tensor_u<double, symmetry_type::symmetric,4,4,4>;
  //for (int i=0;i<4; ++i) {
  //  for (int j=0;j<4; ++j) {
  //    for (int k=0;k<4; ++k) {
  //       cout << sym3_tensor4_t::multiindex(i,j,k) << ", ";
  //    }
  //    cout << endl;
  //  }
  //  cout << endl;
  //}
  
  //cout << "--- Symmetric tensor, rank 4, dimension 3: ---" << endl;
  //using sym_tensor4_t = tensor_u<double, symmetry_type::symmetric,3,3,3,3>;
  //for (int m=0;m<3; ++m) {
  //  for (int i=0;i<3; ++i) {
  //    for (int j=0;j<3; ++j) {
  //      for (int k=0;k<3; ++k) {
  //         cout << sym_tensor4_t::multiindex(m,i,j,k) << ", ";
  //      }
  //      cout << " -1,"<< endl;
  //    }
  //    cout << "-1, -1, -1, -1," << endl
  //         << endl;
  //  }
  //  cout << "-1, -1, -1, -1," << endl
  //       << "-1, -1, -1, -1," << endl
  //       << "-1, -1, -1, -1," << endl
  //       << "-1, -1, -1, -1," << endl
  //       << endl;
  //}

  //cout << "--- Symmetric tensor, rank 4, dimension 2: ---" << endl;
  //using sym_tensor4_t = tensor_u<double, symmetry_type::symmetric,2,2,2,2>;
  //for (int m=0;m<2; ++m) {
  //  for (int i=0;i<2; ++i) {
  //    for (int j=0;j<2; ++j) {
  //      for (int k=0;k<2; ++k) {
  //         cout << sym_tensor4_t::multiindex(m,i,j,k) << ", ";
  //      }
  //      cout << " -1, -1,"<< endl;
  //    }
  //    cout << "-1, -1, -1, -1," << endl
  //         << "-1, -1, -1, -1," << endl
  //         << endl;
  //  }
  //  cout << "-1, -1, -1, -1," << endl
  //       << "-1, -1, -1, -1," << endl
  //       << "-1, -1, -1, -1," << endl
  //       << "-1, -1, -1, -1," << endl
  //       << endl;
  //  cout << "-1, -1, -1, -1," << endl
  //       << "-1, -1, -1, -1," << endl
  //       << "-1, -1, -1, -1," << endl
  //       << "-1, -1, -1, -1," << endl
  //       << endl;
  //}

  //cout << "--- Symmetric tensor, rank 4, dimension 4: ---" << endl;
  //using sym_tensor4_t = tensor_u<double, symmetry_type::symmetric,4,4,4,4>;
  //for (int m=0;m<4; ++m) {
  //  for (int i=0;i<4; ++i) {
  //    for (int j=0;j<4; ++j) {
  //      for (int k=0;k<4; ++k) {
  //         cout << sym_tensor4_t::multiindex(m,i,j,k) << ", ";
  //      }
  //      cout << endl;
  //    }
  //    cout << endl;
  //  }
  //  cout << endl;
  //}
  
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

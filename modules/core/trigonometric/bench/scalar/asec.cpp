//==============================================================================
//         Copyright 2003 - 2013   LASMEA UMR 6602 CNRS/Univ. Clermont II
//         Copyright 2009 - 2013   LRI    UMR 8623 CNRS/Univ Paris Sud XI
//
//          Distributed under the Boost Software License, Version 1.0.
//                 See accompanying file LICENSE.txt or copy at
//                     http://www.boost.org/LICENSE_1_0.txt
//==============================================================================
#include <nt2/trigonometric/include/functions/asec.hpp>

#include <nt2/sdk/bench/benchmark.hpp>
#include <nt2/sdk/bench/timing.hpp>
#include <nt2/include/constants/one.hpp>
#include <nt2/include/constants/ten.hpp>

typedef NT2_SIMD_DEFAULT_EXTENSION  ext_t;
using nt2::tag::asec_;

#define RS(T,V1,V2) (T, (V1) ,(V2))

namespace n1 {
  typedef float T;
  NT2_TIMING(asec_,(RS(T,T(1.01),T(10))))
  NT2_TIMING(asec_,(RS(T,T(-10),T(-1.01))))
}
namespace n2 {
  typedef double T;
  NT2_TIMING(asec_,(RS(T,T(1.01),T(10))))
  NT2_TIMING(asec_,(RS(T,T(-10),T(-1.01))))
}
#undef RS

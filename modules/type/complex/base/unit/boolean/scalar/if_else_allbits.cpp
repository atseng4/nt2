//==============================================================================
//         Copyright 2003 - 2013   LASMEA UMR 6602 CNRS/Univ. Clermont II
//         Copyright 2009 - 2013   LRI    UMR 8623 CNRS/Univ Paris Sud XI
//
//          Distributed under the Boost Software License, Version 1.0.
//                 See accompanying file LICENSE.txt or copy at
//                     http://www.boost.org/LICENSE_1_0.txt
//==============================================================================
#include <nt2/include/functions/if_else_allbits.hpp>

#include <boost/dispatch/functor/meta/call.hpp>
#include <nt2/sdk/functor/meta/call.hpp>
#include <nt2/sdk/unit/tests/relation.hpp>
#include <nt2/sdk/unit/tests/type_expr.hpp>
#include <complex>
#include <nt2/sdk/complex/complex.hpp>
#include <nt2/sdk/unit/module.hpp>
#include <boost/simd/sdk/config.hpp>

#include <nt2/include/constants/false.hpp>
#include <nt2/include/constants/true.hpp>
#include <nt2/include/constants/zero.hpp>
#include <nt2/include/constants/one.hpp>
#include <nt2/include/constants/inf.hpp>
#include <nt2/include/constants/minf.hpp>
#include <nt2/include/constants/nan.hpp>
#include <nt2/include/constants/cnan.hpp>

NT2_TEST_CASE_TPL ( if_else_allbits_real,  NT2_REAL_TYPES)
{
  using nt2::if_else_allbits;
  using nt2::tag::if_else_allbits_;
  typedef typename std::complex<T> cT;
  typedef typename nt2::meta::call<if_else_allbits_(cT, cT)>::type r_t;

  // return type conformity test
  NT2_TEST_TYPE_IS(r_t, cT);

  // specific values tests
#ifndef BOOST_SIMD_NO_INVALIDS
  NT2_TEST_EQUAL(if_else_allbits(nt2::Inf<cT>(), cT(1)), cT(1));
  NT2_TEST_EQUAL(if_else_allbits(nt2::Minf<cT>(), cT(1)),  cT(1));
  NT2_TEST_EQUAL(if_else_allbits(nt2::Nan<cT>(), cT(1)) ,  cT(1));
#endif
  NT2_TEST_EQUAL(if_else_allbits(nt2::False<nt2::logical<T> >(), cT(1)), nt2::Cnan<cT>());
  NT2_TEST_EQUAL(if_else_allbits(nt2::True<nt2::logical<T> >(),cT(1)), cT(1));
  NT2_TEST_EQUAL(if_else_allbits(nt2::Zero<cT>(), cT(1)),  nt2::Cnan<cT>());
  NT2_TEST_EQUAL(if_else_allbits(nt2::One <cT>(), cT(2)),  cT(2));
}

//==============================================================================
//         Copyright 2003 - 2013   LASMEA UMR 6602 CNRS/Univ. Clermont II
//         Copyright 2009 - 2013   LRI    UMR 8623 CNRS/Univ Paris Sud XI
//
//          Distributed under the Boost Software License, Version 1.0.
//                 See accompanying file LICENSE.txt or copy at
//                     http://www.boost.org/LICENSE_1_0.txt
//==============================================================================
#include <nt2/exponential/include/functions/log1p.hpp>

#include <nt2/sdk/functor/meta/call.hpp>
#include <nt2/sdk/unit/module.hpp>
#include <nt2/sdk/unit/tests/relation.hpp>
#include <nt2/sdk/unit/tests/type_expr.hpp>
#include <nt2/sdk/unit/tests/ulp.hpp>
#include <nt2/sdk/unit/module.hpp>
#include <boost/simd/sdk/config.hpp>
#include <nt2/sdk/meta/as_floating.hpp>

#include <nt2/include/constants/eps.hpp>
#include <nt2/include/constants/mone.hpp>
#include <nt2/include/constants/one.hpp>
#include <nt2/include/constants/smallestposval.hpp>
#include <nt2/include/constants/zero.hpp>
#include <nt2/include/constants/minf.hpp>
#include <nt2/include/constants/nan.hpp>
#include <nt2/include/constants/log_2.hpp>

NT2_TEST_CASE_TPL ( log1p_real,  NT2_REAL_TYPES)
{
  using nt2::log1p;
  using nt2::tag::log1p_;

  typedef typename nt2::meta::call<log1p_(T)>::type r_t;
  typedef T wished_r_t;

  // return type conformity test
  NT2_TEST_TYPE_IS(r_t, wished_r_t);

  // specific values tests
#ifndef BOOST_SIMD_NO_INVALIDS
  NT2_TEST_ULP_EQUAL(log1p(nt2::Nan<T>()), nt2::Nan<r_t>(), 0.5);
  NT2_TEST_ULP_EQUAL(log1p(nt2::Minf<T>()), nt2::Nan<r_t>(), 0.5);
  NT2_TEST_ULP_EQUAL(log1p(nt2::Mone<T>()), nt2::Minf<r_t>(), 0.5);
#endif
  NT2_TEST_ULP_EQUAL(log1p(nt2::Eps<T>()), nt2::Eps<r_t>(), 0.5);
  NT2_TEST_ULP_EQUAL(log1p(nt2::Eps<T>()), nt2::Eps<r_t>(), 0.5);
  NT2_TEST_ULP_EQUAL(log1p(nt2::One<T>()), nt2::Log_2<r_t>(), 0.5);
  NT2_TEST_ULP_EQUAL(log1p(nt2::Zero<T>()), nt2::Zero<r_t>(), 0.5);
  NT2_TEST_ULP_EQUAL(log1p(nt2::Smallestposval<T>()), nt2::Smallestposval<T>(), 0.5);
}

NT2_TEST_CASE_TPL ( log1p_unsigned_int,  NT2_UNSIGNED_TYPES)
{
  using nt2::log1p;
  using nt2::tag::log1p_;

  typedef typename nt2::meta::call<log1p_(T)>::type r_t;
  typedef typename nt2::meta::as_floating<T>::type wished_r_t;

  // return type conformity test
  NT2_TEST_TYPE_IS(r_t, wished_r_t);

  // specific values tests
  NT2_TEST_ULP_EQUAL(log1p(nt2::One<T>()), nt2::Log_2<r_t>(), 0.5);
  NT2_TEST_ULP_EQUAL(log1p(nt2::Zero<T>()), nt2::Zero<r_t>(), 0.5);
}

NT2_TEST_CASE_TPL ( log1p_signed_int,  NT2_INTEGRAL_SIGNED_TYPES)
{
  using nt2::log1p;
  using nt2::tag::log1p_;

  typedef typename nt2::meta::call<log1p_(T)>::type r_t;
  typedef typename nt2::meta::as_floating<T>::type wished_r_t;

  // return type conformity test
  NT2_TEST_TYPE_IS(r_t, wished_r_t);

  // specific values tests
  NT2_TEST_ULP_EQUAL(log1p(nt2::Mone<T>()), nt2::Minf<r_t>(), 0.5);
  NT2_TEST_ULP_EQUAL(log1p(nt2::One<T>()), nt2::Log_2<r_t>(), 0.5);
  NT2_TEST_ULP_EQUAL(log1p(nt2::Zero<T>()), nt2::Zero<r_t>(), 0.5);
}

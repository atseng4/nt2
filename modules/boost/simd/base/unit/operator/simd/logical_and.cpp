//==============================================================================
//         Copyright 2003 - 2012   LASMEA UMR 6602 CNRS/Univ. Clermont II
//         Copyright 2009 - 2012   LRI    UMR 8623 CNRS/Univ Paris Sud XI
//
//          Distributed under the Boost Software License, Version 1.0.
//                 See accompanying file LICENSE.txt or copy at
//                     http://www.boost.org/LICENSE_1_0.txt
//==============================================================================
#include <boost/simd/operator/include/functions/logical_and.hpp>
#include <boost/simd/sdk/simd/native.hpp>
#include <boost/simd/sdk/simd/logical.hpp>

#include <boost/dispatch/functor/meta/call.hpp>
#include <boost/simd/sdk/simd/native.hpp>
#include <nt2/sdk/unit/tests/relation.hpp>
#include <nt2/sdk/unit/tests/type_expr.hpp>
#include <nt2/sdk/unit/module.hpp>
#include <boost/simd/sdk/config.hpp>
#include <boost/simd/sdk/simd/io.hpp>

#include <boost/simd/include/constants/one.hpp>
#include <boost/simd/include/constants/half.hpp>
#include <boost/simd/include/constants/mone.hpp>
#include <boost/simd/include/constants/zero.hpp>
#include <boost/simd/include/constants/inf.hpp>
#include <boost/simd/include/constants/minf.hpp>
#include <boost/simd/include/constants/nan.hpp>
#include <boost/simd/include/constants/false.hpp>
#include <boost/simd/include/constants/true.hpp>

NT2_TEST_CASE_TPL ( logical_and_integer,  BOOST_SIMD_SIMD_INTEGRAL_TYPES)
{
  using boost::simd::logical_and;
  using boost::simd::tag::logical_and_;
  using boost::simd::native;
  typedef BOOST_SIMD_DEFAULT_EXTENSION  ext_t;
  typedef native<T,ext_t>                  vT;
  typedef native< boost::simd::logical<T>, ext_t> vlT;
  typedef typename boost::dispatch::meta::call<logical_and_(vT,vT)>::type r_t;

  // specific values tests
  NT2_TEST_EQUAL(logical_and(boost::simd::True<vlT>(), boost::simd::True<vlT>()), boost::simd::True<r_t>());
  NT2_TEST_EQUAL(logical_and(boost::simd::True<vlT>(),boost::simd::False<vlT>()), boost::simd::False<r_t>());
  NT2_TEST_EQUAL(logical_and(boost::simd::False<vlT>(), boost::simd::False<vlT>()), boost::simd::False<r_t>());
  NT2_TEST_EQUAL(logical_and(boost::simd::One<vT>(), boost::simd::One<vT>()), boost::simd::True<r_t>());
} // end of test for integer_

NT2_TEST_CASE_TPL ( logical_and_real,  BOOST_SIMD_SIMD_REAL_TYPES)
{
  using boost::simd::logical_and;
  using boost::simd::tag::logical_and_;
  using boost::simd::native;
  typedef BOOST_SIMD_DEFAULT_EXTENSION  ext_t;
  typedef native<T,ext_t>                  vT;
  typedef native< boost::simd::logical<T>, ext_t> vlT;
  typedef typename boost::dispatch::meta::call<logical_and_(vT,vT)>::type r_t;

  // specific values tests
  NT2_TEST_EQUAL(logical_and(boost::simd::Inf<vT>(), boost::simd::Inf<vT>()), boost::simd::True<r_t>());
  NT2_TEST_EQUAL(logical_and(boost::simd::Minf<vT>(), boost::simd::Minf<vT>()), boost::simd::True<r_t>());
  NT2_TEST_EQUAL(logical_and(boost::simd::Nan<vT>(), boost::simd::Nan<vT>()), boost::simd::True<r_t>());
  NT2_TEST_EQUAL(logical_and(boost::simd::One<vT>(),boost::simd::Zero<vT>()), boost::simd::False<r_t>());
  NT2_TEST_EQUAL(logical_and(boost::simd::Zero<vT>(), boost::simd::Zero<vT>()), boost::simd::False<r_t>());
  NT2_TEST_EQUAL(logical_and(boost::simd::Half<vT>(), boost::simd::Half<vT>()), boost::simd::True<r_t>());
} // end of test for floating_

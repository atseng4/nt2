//////////////////////////////////////////////////////////////////////////////
///   Copyright 2003 and onward LASMEA UMR 6602 CNRS/U.B.P Clermont-Ferrand
///   Copyright 2009 and onward LRI    UMR 8623 CNRS/Univ Paris Sud XI
///
///          Distributed under the Boost Software License, Version 1.0
///                 See accompanying file LICENSE.txt or copy at
///                     http://www.boost.org/LICENSE_1_0.txt
//////////////////////////////////////////////////////////////////////////////
#define NT2_UNIT_MODULE "nt2 arithmetic toolbox - toint/simd Mode"

//////////////////////////////////////////////////////////////////////////////
// unit test behavior of arithmetic components in simd mode
//////////////////////////////////////////////////////////////////////////////
/// created by jt the 04/12/2010
/// 
#include <nt2/toolbox/arithmetic/include/functions/toint.hpp>
#include <nt2/include/functions/ulpdist.hpp>
#include <boost/type_traits/is_same.hpp>
#include <nt2/sdk/functor/meta/call.hpp>
#include <nt2/sdk/meta/as_integer.hpp>
#include <nt2/sdk/meta/as_floating.hpp>
#include <nt2/sdk/meta/as_signed.hpp>
#include <nt2/sdk/meta/upgrade.hpp>
#include <nt2/sdk/meta/downgrade.hpp>
#include <nt2/sdk/meta/scalar_of.hpp>
#include <boost/dispatch/meta/as_floating.hpp>
#include <boost/type_traits/common_type.hpp>
#include <nt2/sdk/unit/tests.hpp>
#include <nt2/sdk/unit/module.hpp>

#include <nt2/toolbox/constant/constant.hpp>
#include <nt2/sdk/meta/cardinal_of.hpp>
#include <nt2/include/functions/splat.hpp>

#include <nt2/include/functions/load.hpp>


NT2_TEST_CASE_TPL ( toint_real__1_0,  NT2_SIMD_REAL_TYPES)
{
  using nt2::toint;
  using nt2::tag::toint_;
  using nt2::load; 
  using boost::simd::native;
  using nt2::meta::cardinal_of;
  typedef NT2_SIMD_DEFAULT_EXTENSION  ext_t;
  typedef typename nt2::meta::upgrade<T>::type   u_t;
  typedef native<T,ext_t>                        n_t;
  typedef n_t                                     vT;
  typedef typename nt2::meta::as_integer<T>::type iT;
  typedef native<iT,ext_t>                       ivT;
  typedef typename nt2::meta::strip<typename nt2::meta::call<toint_(vT)>::type>::type r_t;
  typedef typename nt2::meta::strip<typename nt2::meta::call<toint_(T)>::type>::type sr_t;
  typedef typename nt2::meta::scalar_of<r_t>::type ssr_t;
  double ulpd;
  ulpd=0.0;


  // specific values tests
  NT2_TEST_EQUAL(toint(nt2::Inf<vT>())[0], nt2::Inf<sr_t>());
  NT2_TEST_EQUAL(toint(nt2::Minf<vT>())[0], nt2::Minf<sr_t>());
  NT2_TEST_EQUAL(toint(nt2::Mone<vT>())[0], nt2::Mone<sr_t>());
  NT2_TEST_EQUAL(toint(nt2::Nan<vT>())[0], nt2::Zero<sr_t>());
  NT2_TEST_EQUAL(toint(nt2::One<vT>())[0], nt2::One<sr_t>());
  NT2_TEST_EQUAL(toint(nt2::Zero<vT>())[0], nt2::Zero<sr_t>());
} // end of test for floating_
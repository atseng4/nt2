//////////////////////////////////////////////////////////////////////////////
///   Copyright 2003 and onward LASMEA UMR 6602 CNRS/U.B.P Clermont-Ferrand
///   Copyright 2009 and onward LRI    UMR 8623 CNRS/Univ Paris Sud XI
///
///          Distributed under the Boost Software License, Version 1.0
///                 See accompanying file LICENSE.txt or copy at
///                     http://www.boost.org/LICENSE_1_0.txt
//////////////////////////////////////////////////////////////////////////////
#define NT2_UNIT_MODULE "nt2 bitwise toolbox - negif/simd Mode"

//////////////////////////////////////////////////////////////////////////////
// unit test behavior of bitwise components in simd mode
//////////////////////////////////////////////////////////////////////////////
/// created  by jt the 18/02/2011
/// 
#include <nt2/toolbox/bitwise/include/negif.hpp>
#include <nt2/include/functions/ulpdist.hpp>
#include <boost/type_traits/is_same.hpp>
#include <nt2/sdk/functor/meta/call.hpp>
#include <nt2/sdk/unit/tests.hpp>
#include <nt2/sdk/unit/module.hpp>
#include <nt2/sdk/memory/buffer.hpp>
#include <nt2/include/constants/real.hpp>
#include <nt2/include/constants/infinites.hpp>
#include <nt2/sdk/memory/is_aligned.hpp>
#include <nt2/sdk/memory/aligned_type.hpp>
#include <nt2/include/functions/load.hpp>


NT2_TEST_CASE_TPL ( negif_real__2_0,  BOOST_SIMD_SIMD_REAL_TYPES)
{
  using nt2::negif;
  using nt2::tag::negif_;
  using nt2::load; 
  using boost::simd::native;
  using nt2::meta::cardinal_of;
  typedef BOOST_SIMD_DEFAULT_EXTENSION  ext_t;
  typedef typename nt2::meta::upgrade<T>::type   u_t;
  typedef native<T,ext_t>                        n_t;
  typedef n_t                                     vT;
  typedef typename nt2::meta::as_integer<T>::type iT;
  typedef native<iT,ext_t>                       ivT;
  typedef typename nt2::meta::call<negif_(vT,vT)>::type r_t;
  typedef typename nt2::meta::call<negif_(T,T)>::type sr_t;
  typedef typename nt2::meta::scalar_of<r_t>::type ssr_t;
  double ulpd;
  ulpd=0.0;


  // specific values tests
  NT2_TEST_EQUAL(negif(nt2::splat<vT>(0),nt2::splat<vT>(1))[0], 1);
  NT2_TEST_EQUAL(negif(nt2::splat<vT>(1),nt2::splat<vT>(1))[0], -1);
  NT2_TEST_EQUAL(negif(nt2::Inf<vT>(),nt2::splat<vT>(1))[0], -1);
  NT2_TEST_EQUAL(negif(nt2::Minf<vT>(),nt2::splat<vT>(1))[0], -1);
  NT2_TEST_EQUAL(negif(nt2::Nan<vT>(),nt2::splat<vT>(1))[0], -1);
  NT2_TEST_EQUAL(negif(nt2::Zero<vT>(),nt2::splat<vT>(1))[0], 1);
} // end of test for real_

NT2_TEST_CASE_TPL ( negif_signed_int__2_0,  BOOST_SIMD_SIMD_INTEGRAL_SIGNED_TYPES)
{
  using nt2::negif;
  using nt2::tag::negif_;
  using nt2::load; 
  using boost::simd::native;
  using nt2::meta::cardinal_of;
  typedef BOOST_SIMD_DEFAULT_EXTENSION  ext_t;
  typedef typename nt2::meta::upgrade<T>::type   u_t;
  typedef native<T,ext_t>                        n_t;
  typedef n_t                                     vT;
  typedef typename nt2::meta::as_integer<T>::type iT;
  typedef native<iT,ext_t>                       ivT;
  typedef typename nt2::meta::call<negif_(vT,vT)>::type r_t;
  typedef typename nt2::meta::call<negif_(T,T)>::type sr_t;
  typedef typename nt2::meta::scalar_of<r_t>::type ssr_t;
  double ulpd;
  ulpd=0.0;


  // specific values tests
  NT2_TEST_EQUAL(negif(nt2::splat<vT>(0),nt2::splat<vT>(1))[0], 1);
  NT2_TEST_EQUAL(negif(nt2::splat<vT>(1),nt2::splat<vT>(1))[0], sr_t(-1));
  NT2_TEST_EQUAL(negif(nt2::Zero<vT>(), nt2::Zero<vT>())[0], nt2::Zero<sr_t>());
} // end of test for signed_int_

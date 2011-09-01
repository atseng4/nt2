//////////////////////////////////////////////////////////////////////////////
///   Copyright 2003 and onward LASMEA UMR 6602 CNRS/U.B.P Clermont-Ferrand
///   Copyright 2009 and onward LRI    UMR 8623 CNRS/Univ Paris Sud XI
///
///          Distributed under the Boost Software License, Version 1.0
///                 See accompanying file LICENSE.txt or copy at
///                     http://www.boost.org/LICENSE_1_0.txt
//////////////////////////////////////////////////////////////////////////////
#define NT2_UNIT_MODULE "nt2 boost.simd.reduction toolbox - hmsb/simd Mode"

//////////////////////////////////////////////////////////////////////////////
// cover test behavior of boost.simd.reduction components in simd mode
//////////////////////////////////////////////////////////////////////////////
/// created  by jt the 24/02/2011
/// 
#include <boost/simd/toolbox/reduction/include/functions/hmsb.hpp>
#include <boost/simd/include/functions/ulpdist.hpp>
#include <boost/simd/include/functions/max.hpp>
#include <boost/simd/include/functions/bits.hpp>
#include <boost/simd/include/functions/shri.hpp>

#include <boost/simd/include/functions/bits.hpp>

#include <boost/type_traits/is_same.hpp>
#include <boost/dispatch/functor/meta/call.hpp>
#include <nt2/sdk/unit/tests.hpp>
#include <nt2/sdk/unit/module.hpp>
#include <boost/simd/sdk/memory/buffer.hpp>
#include <boost/simd/toolbox/constant/constant.hpp>
#include <boost/simd/sdk/memory/is_aligned.hpp>
#include <boost/simd/sdk/memory/aligned_type.hpp>
#include <boost/simd/include/functions/load.hpp>
#include <boost/simd/toolbox/constant/constant.hpp>


NT2_TEST_CASE_TPL ( hmsb_real__1_0,  BOOST_SIMD_SIMD_REAL_TYPES)
{
  using boost::simd::hmsb;
  using boost::simd::tag::hmsb_;
  using boost::simd::load; 
  using boost::simd::native;
  using boost::simd::meta::cardinal_of;
  typedef BOOST_SIMD_DEFAULT_EXTENSION  ext_t;
  typedef typename boost::dispatch::meta::upgrade<T>::type   u_t;
  typedef native<T,ext_t>                        n_t;
  typedef n_t                                     vT;
  typedef typename boost::dispatch::meta::as_integer<T>::type iT;
  typedef native<iT,ext_t>                       ivT;
  typedef typename boost::dispatch::meta::call<hmsb_(vT)>::type r_t;
  typedef typename boost::simd::meta::scalar_of<r_t>::type sr_t;
  typedef typename boost::simd::meta::scalar_of<r_t>::type ssr_t;
  double ulpd;
  ulpd=0.0;

  // random verifications
  static const nt2::uint32_t NR = NT2_NB_RANDOM_TEST;
  {
    NT2_CREATE_BUF(tab_a0,T, NR, boost::simd::Valmin<T>(), boost::simd::Valmax<T>());
    double ulp0, ulpd ; ulpd=ulp0=0.0;
    for(nt2::uint32_t j = 0; j < NR/cardinal_of<n_t>::value; j++)
      {
        vT a0 = load<vT>(&tab_a0[0],j);
        r_t v = hmsb(a0);
        int z = 0;
        int N = cardinal_of<n_t>::value;
        for(int i = 0; i< N; ++i)
        {
          z |= boost::simd::bits(a0[i]) >> (sizeof(iT)*CHAR_BIT - 1) << (N-i-1);        }
        NT2_TEST_EQUAL( v,z);
      }
    
  }
} // end of test for real_

NT2_TEST_CASE_TPL ( hmsb_uint64_t_1_0,  (nt2::uint64_t))
{
  using boost::simd::hmsb;
  using boost::simd::tag::hmsb_;
  using boost::simd::load; 
  using boost::simd::native;
  using boost::simd::meta::cardinal_of;
  typedef BOOST_SIMD_DEFAULT_EXTENSION  ext_t;
  typedef typename boost::dispatch::meta::upgrade<T>::type   u_t;
  typedef native<T,ext_t>                        n_t;
  typedef n_t                                     vT;
  typedef typename boost::dispatch::meta::as_integer<T>::type iT;
  typedef native<iT,ext_t>                       ivT;
  typedef typename boost::dispatch::meta::call<hmsb_(vT)>::type r_t;
  typedef typename boost::simd::meta::scalar_of<r_t>::type sr_t;
  typedef typename boost::simd::meta::scalar_of<r_t>::type ssr_t;
  double ulpd;
  ulpd=0.0;

  // random verifications
  static const nt2::uint32_t NR = NT2_NB_RANDOM_TEST;
  {
    NT2_CREATE_BUF(tab_a0,T, NR, boost::simd::Valmin<T>(), boost::simd::Valmax<T>());
    double ulp0, ulpd ; ulpd=ulp0=0.0;
    for(nt2::uint32_t j = 0; j < NR/cardinal_of<n_t>::value; j++)
      {
        vT a0 = load<vT>(&tab_a0[0],j);
        r_t v = hmsb(a0);
        int z = 0;
        int N = cardinal_of<n_t>::value;
        for(int i = 0; i< N; ++i)
        {
          z |= boost::simd::bits(a0[i]) >> (sizeof(iT)*CHAR_BIT - 1) << (N-i-1);        }
        NT2_TEST_EQUAL( v,z);
      }
    
  }
} // end of test for uint64_t

NT2_TEST_CASE_TPL ( hmsb_int64_t_1_0,  (nt2::int64_t))
{
  using boost::simd::hmsb;
  using boost::simd::tag::hmsb_;
  using boost::simd::load; 
  using boost::simd::native;
  using boost::simd::meta::cardinal_of;
  typedef BOOST_SIMD_DEFAULT_EXTENSION  ext_t;
  typedef typename boost::dispatch::meta::upgrade<T>::type   u_t;
  typedef native<T,ext_t>                        n_t;
  typedef n_t                                     vT;
  typedef typename boost::dispatch::meta::as_integer<T>::type iT;
  typedef native<iT,ext_t>                       ivT;
  typedef typename boost::dispatch::meta::call<hmsb_(vT)>::type r_t;
  typedef typename boost::simd::meta::scalar_of<r_t>::type sr_t;
  typedef typename boost::simd::meta::scalar_of<r_t>::type ssr_t;
  double ulpd;
  ulpd=0.0;

  // random verifications
  static const nt2::uint32_t NR = NT2_NB_RANDOM_TEST;
  {
    NT2_CREATE_BUF(tab_a0,T, NR, boost::simd::Valmin<T>(), boost::simd::Valmax<T>());
    double ulp0, ulpd ; ulpd=ulp0=0.0;
    for(nt2::uint32_t j = 0; j < NR/cardinal_of<n_t>::value; j++)
      {
        vT a0 = load<vT>(&tab_a0[0],j);
        r_t v = hmsb(a0);
        int z = 0;
        int N = cardinal_of<n_t>::value;
        for(int i = 0; i< N; ++i)
        {
          z |= boost::simd::bits(a0[i]) >> (sizeof(iT)*CHAR_BIT - 1) << (N-i-1);        }
        NT2_TEST_EQUAL( v,z);
      }
    
  }
} // end of test for int64_t

NT2_TEST_CASE_TPL ( hmsb_uint32_t_1_0,  (nt2::uint32_t))
{
  using boost::simd::hmsb;
  using boost::simd::tag::hmsb_;
  using boost::simd::load; 
  using boost::simd::native;
  using boost::simd::meta::cardinal_of;
  typedef BOOST_SIMD_DEFAULT_EXTENSION  ext_t;
  typedef typename boost::dispatch::meta::upgrade<T>::type   u_t;
  typedef native<T,ext_t>                        n_t;
  typedef n_t                                     vT;
  typedef typename boost::dispatch::meta::as_integer<T>::type iT;
  typedef native<iT,ext_t>                       ivT;
  typedef typename boost::dispatch::meta::call<hmsb_(vT)>::type r_t;
  typedef typename boost::simd::meta::scalar_of<r_t>::type sr_t;
  typedef typename boost::simd::meta::scalar_of<r_t>::type ssr_t;
  double ulpd;
  ulpd=0.0;

  // random verifications
  static const nt2::uint32_t NR = NT2_NB_RANDOM_TEST;
  {
    NT2_CREATE_BUF(tab_a0,T, NR, boost::simd::Valmin<T>(), boost::simd::Valmax<T>());
    double ulp0, ulpd ; ulpd=ulp0=0.0;
    for(nt2::uint32_t j = 0; j < NR/cardinal_of<n_t>::value; j++)
      {
        vT a0 = load<vT>(&tab_a0[0],j);
        r_t v = hmsb(a0);
        int z = 0;
        int N = cardinal_of<n_t>::value;
        for(int i = 0; i< N; ++i)
        {
          z |= boost::simd::bits(a0[i]) >> (sizeof(iT)*CHAR_BIT - 1) << (N-i-1);        }
        NT2_TEST_EQUAL( v,z);
      }
    
  }
} // end of test for uint32_t

NT2_TEST_CASE_TPL ( hmsb_int32_t_1_0,  (nt2::int32_t))
{
  using boost::simd::hmsb;
  using boost::simd::tag::hmsb_;
  using boost::simd::load; 
  using boost::simd::native;
  using boost::simd::meta::cardinal_of;
  typedef BOOST_SIMD_DEFAULT_EXTENSION  ext_t;
  typedef typename boost::dispatch::meta::upgrade<T>::type   u_t;
  typedef native<T,ext_t>                        n_t;
  typedef n_t                                     vT;
  typedef typename boost::dispatch::meta::as_integer<T>::type iT;
  typedef native<iT,ext_t>                       ivT;
  typedef typename boost::dispatch::meta::call<hmsb_(vT)>::type r_t;
  typedef typename boost::simd::meta::scalar_of<r_t>::type sr_t;
  typedef typename boost::simd::meta::scalar_of<r_t>::type ssr_t;
  double ulpd;
  ulpd=0.0;

  // random verifications
  static const nt2::uint32_t NR = NT2_NB_RANDOM_TEST;
  {
    NT2_CREATE_BUF(tab_a0,T, NR, boost::simd::Valmin<T>(), boost::simd::Valmax<T>());
    double ulp0, ulpd ; ulpd=ulp0=0.0;
    for(nt2::uint32_t j = 0; j < NR/cardinal_of<n_t>::value; j++)
      {
        vT a0 = load<vT>(&tab_a0[0],j);
        r_t v = hmsb(a0);
        int z = 0;
        int N = cardinal_of<n_t>::value;
        for(int i = 0; i< N; ++i)
        {
          z |= boost::simd::bits(a0[i]) >> (sizeof(iT)*CHAR_BIT - 1) << (N-i-1);        }
        NT2_TEST_EQUAL( v,z);
      }
    
  }
} // end of test for int32_t

NT2_TEST_CASE_TPL ( hmsb_uint8_t_1_0,  (nt2::uint8_t))
{
  using boost::simd::hmsb;
  using boost::simd::tag::hmsb_;
  using boost::simd::load; 
  using boost::simd::native;
  using boost::simd::meta::cardinal_of;
  typedef BOOST_SIMD_DEFAULT_EXTENSION  ext_t;
  typedef typename boost::dispatch::meta::upgrade<T>::type   u_t;
  typedef native<T,ext_t>                        n_t;
  typedef n_t                                     vT;
  typedef typename boost::dispatch::meta::as_integer<T>::type iT;
  typedef native<iT,ext_t>                       ivT;
  typedef typename boost::dispatch::meta::call<hmsb_(vT)>::type r_t;
  typedef typename boost::simd::meta::scalar_of<r_t>::type sr_t;
  typedef typename boost::simd::meta::scalar_of<r_t>::type ssr_t;
  double ulpd;
  ulpd=0.0;

  // random verifications
  static const nt2::uint32_t NR = NT2_NB_RANDOM_TEST;
  {
    NT2_CREATE_BUF(tab_a0,T, NR, boost::simd::Valmin<T>(), boost::simd::Valmax<T>());
    double ulp0, ulpd ; ulpd=ulp0=0.0;
    for(nt2::uint32_t j = 0; j < NR/cardinal_of<n_t>::value; j++)
      {
        vT a0 = load<vT>(&tab_a0[0],j);
        r_t v = hmsb(a0);
        int z = 0;
        int N = cardinal_of<n_t>::value;
        for(int i = 0; i< N; ++i)
        {
          z |= boost::simd::bits(a0[i]) >> (sizeof(iT)*CHAR_BIT - 1) << (N-i-1);        }
        NT2_TEST_EQUAL( v,z);
      }
    
  }
} // end of test for uint8_t

NT2_TEST_CASE_TPL ( hmsb_int8_t_1_0,  (nt2::int8_t))
{
  using boost::simd::hmsb;
  using boost::simd::tag::hmsb_;
  using boost::simd::load; 
  using boost::simd::native;
  using boost::simd::meta::cardinal_of;
  typedef BOOST_SIMD_DEFAULT_EXTENSION  ext_t;
  typedef typename boost::dispatch::meta::upgrade<T>::type   u_t;
  typedef native<T,ext_t>                        n_t;
  typedef n_t                                     vT;
  typedef typename boost::dispatch::meta::as_integer<T>::type iT;
  typedef native<iT,ext_t>                       ivT;
  typedef typename boost::dispatch::meta::call<hmsb_(vT)>::type r_t;
  typedef typename boost::simd::meta::scalar_of<r_t>::type sr_t;
  typedef typename boost::simd::meta::scalar_of<r_t>::type ssr_t;
  double ulpd;
  ulpd=0.0;

  // random verifications
  static const nt2::uint32_t NR = NT2_NB_RANDOM_TEST;
  {
    NT2_CREATE_BUF(tab_a0,T, NR, boost::simd::Valmin<T>(), boost::simd::Valmax<T>());
    double ulp0, ulpd ; ulpd=ulp0=0.0;
    for(nt2::uint32_t j = 0; j < NR/cardinal_of<n_t>::value; j++)
      {
        vT a0 = load<vT>(&tab_a0[0],j);
        r_t v = hmsb(a0);
        int z = 0;
        int N = cardinal_of<n_t>::value;
        for(int i = 0; i< N; ++i)
        {
          z |= boost::simd::bits(a0[i]) >> (sizeof(iT)*CHAR_BIT - 1) << (N-i-1);        }
        NT2_TEST_EQUAL( v,z);
      }
    
  }
} // end of test for int8_t

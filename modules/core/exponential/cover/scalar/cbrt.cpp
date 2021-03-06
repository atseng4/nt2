//==============================================================================
//         Copyright 2003 - 2012   LASMEA UMR 6602 CNRS/Univ. Clermont II
//         Copyright 2009 - 2012   LRI    UMR 8623 CNRS/Univ Paris Sud XI
//
//          Distributed under the Boost Software License, Version 1.0.
//                 See accompanying file LICENSE.txt or copy at
//                     http://www.boost.org/LICENSE_1_0.txt
//==============================================================================
#define NT2_UNIT_MODULE "nt2 exponential toolbox - cbrt/scalar Mode"

//////////////////////////////////////////////////////////////////////////////
// cover test behavior of exponential components in scalar mode
//////////////////////////////////////////////////////////////////////////////
/// created by jt the 08/12/2010
///
#include <nt2/exponential/include/functions/cbrt.hpp>
#include <nt2/include/functions/sqr.hpp>
extern "C" { long double cephes_cbrtl(long double);}

#include <nt2/sdk/unit/tests/ulp.hpp>
#include <nt2/sdk/unit/tests/cover.hpp>
#include <nt2/sdk/unit/module.hpp>

#include <nt2/constant/constant.hpp>


NT2_TEST_CASE_TPL ( cbrt_real,  NT2_REAL_TYPES)
{

  using nt2::cbrt;
  using nt2::tag::cbrt_;

  // random verifications
  static const nt2::uint32_t NR = NT2_NB_RANDOM_TEST;
  {
    NT2_CREATE_BUF(tab_a0,T, NR, T(-10), T(10));
    T a0;
    for(nt2::uint32_t j =0; j < NR; ++j )
      {
        std::cout << "for param "
                  << "  a0 = "<< (a0 = tab_a0[j])
                  << std::endl;
        NT2_TEST_ULP_EQUAL( nt2::cbrt(a0),T(::cephes_cbrtl(a0)),2);
        NT2_TEST_ULP_EQUAL( nt2::sqr(nt2::cbrt(a0))*nt2::cbrt(a0),a0,3);
     }
   }
} // end of test for floating_


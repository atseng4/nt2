//////////////////////////////////////////////////////////////////////////////
///   Copyright 2003 and onward LASMEA UMR 6602 CNRS/U.B.P Clermont-Ferrand
///   Copyright 2009 and onward LRI    UMR 8623 CNRS/Univ Paris Sud XI
///
///          Distributed under the Boost Software License, Version 1.0
///                 See accompanying file LICENSE.txt or copy at
///                     http://www.boost.org/LICENSE_1_0.txt
//////////////////////////////////////////////////////////////////////////////
#define NT2_UNIT_MODULE "nt2 ieee toolbox - unit/scalar Mode"
#include <nt2/sdk/functor/meta/call.hpp>
#include <boost/type_traits/is_same.hpp>
#include <nt2/toolbox/ieee/include/nextpow2.hpp>
#include <nt2/sdk/unit/tests.hpp>
#include <nt2/sdk/unit/module.hpp>
#include <nt2/sdk/constant/real.hpp>
#include <nt2/sdk/constant/eps_related.hpp>
#include <nt2/sdk/meta/as_real.hpp>
#include <nt2/sdk/meta/as_integer.hpp>

//////////////////////////////////////////////////////////////////////////////
// Test behavior of arithmetic components using NT2_TEST_CASE
//////////////////////////////////////////////////////////////////////////////

NT2_TEST_CASE_TPL ( integral_nextpow2,   NT2_TYPES        
                  )
{
  using nt2::nextpow2;
  using nt2::tag::nextpow2_;
  typedef typename boost::result_of<nt2::meta::floating(T)>::type t1_t; 
  typedef typename nt2::meta::as_integer<t1_t, signed > ::type  t2_t; 
  NT2_TEST( (boost::is_same < typename nt2::meta::call<nextpow2_(T)>::type
           , t2_t
              >::value)
           );
  NT2_TEST_EQUAL(  nextpow2( T(42) ), 6 );
  NT2_TEST_EQUAL(  nextpow2( T(32) ), 5 );
  
    

}

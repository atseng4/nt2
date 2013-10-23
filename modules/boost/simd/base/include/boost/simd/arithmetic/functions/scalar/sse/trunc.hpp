//==============================================================================
//         Copyright 2003 - 2011 LASMEA UMR 6602 CNRS/Univ. Clermont II
//         Copyright 2009 - 2011 LRI    UMR 8623 CNRS/Univ Paris Sud XI
//
//          Distributed under the Boost Software License, Version 1.0.
//                 See accompanying file LICENSE.txt or copy at
//                     http://www.boost.org/LICENSE_1_0.txt
//==============================================================================
#ifndef BOOST_SIMD_ARITHMETIC_FUNCTIONS_SCALAR_SSE_TRUNC_HPP_INCLUDED
#define BOOST_SIMD_ARITHMETIC_FUNCTIONS_SCALAR_SSE_TRUNC_HPP_INCLUDED

#if defined(BOOST_SIMD_HAS_SSE2_SUPPORT)
#include <boost/simd/arithmetic/functions/trunc.hpp>
#include <boost/simd/include/functions/scalar/abs.hpp>
#include <boost/simd/include/functions/splat.hpp>
#include <boost/simd/include/constants/maxflint.hpp>
#include <boost/simd/sdk/simd/native.hpp>

namespace boost { namespace simd { namespace ext
{
  BOOST_SIMD_FUNCTOR_IMPLEMENTATION( boost::simd::tag::trunc_
                                   , boost::simd::tag::sse_
                                   , (A0)
                                   , (scalar_< single_<A0> >)
                                   )
  {

    typedef A0 result_type;
    BOOST_SIMD_FUNCTOR_CALL(1)
    {
      typedef native<result_type, tag::sse_> v_type;
      if (!(abs(a0) <=  Maxflint<result_type>()))
        return a0;
      return result_type(_mm_cvttss_si32(splat<v_type>(a0)));
    }
  };


} } }


#endif
#endif

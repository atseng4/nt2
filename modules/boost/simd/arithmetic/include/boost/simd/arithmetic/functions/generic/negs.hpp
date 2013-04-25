//==============================================================================
//         Copyright 2003 - 2011 LASMEA UMR 6602 CNRS/Univ. Clermont II
//         Copyright 2009 - 2011 LRI    UMR 8623 CNRS/Univ Paris Sud XI
//
//          Distributed under the Boost Software License, Version 1.0.
//                 See accompanying file LICENSE.txt or copy at
//                     http://www.boost.org/LICENSE_1_0.txt
//==============================================================================
#ifndef BOOST_SIMD_ARITHMETIC_FUNCTIONS_GENERIC_NEGS_HPP_INCLUDED
#define BOOST_SIMD_ARITHMETIC_FUNCTIONS_GENERIC_NEGS_HPP_INCLUDED
#include <boost/simd/arithmetic/functions/negs.hpp>
#include <boost/dispatch/meta/as_floating.hpp>
#include <boost/simd/include/functions/simd/if_else.hpp>
#include <boost/simd/include/functions/simd/unary_minus.hpp>

namespace boost { namespace simd { namespace ext
{
  BOOST_SIMD_FUNCTOR_IMPLEMENTATION( boost::simd::tag::negs_, tag::cpu_
                            , (A0)
                            , ((generic_<signed_<A0> >))
                            )
  {
    typedef A0 result_type;
    BOOST_SIMD_FUNCTOR_CALL(1)
    {
      return boost::simd::if_else(eq(a0, boost::simd::Valmin<A0>()),
                                  boost::simd::Valmax<A0>(),
                                  boost::simd::unary_minus(a0));
    }
  };

  BOOST_SIMD_FUNCTOR_IMPLEMENTATION( boost::simd::tag::negs_, tag::cpu_
                            , (A0)
                            , ((generic_<floating_<A0> >))
                            )
  {
    typedef A0 result_type;
    BOOST_SIMD_FUNCTOR_CALL(1) { return boost::simd::unary_minus(a0); }
  };
} } }


#endif

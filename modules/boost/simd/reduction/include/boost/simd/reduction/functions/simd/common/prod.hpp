//==============================================================================
//         Copyright 2003 - 2011 LASMEA UMR 6602 CNRS/Univ. Clermont II
//         Copyright 2009 - 2011 LRI    UMR 8623 CNRS/Univ Paris Sud XI
//
//          Distributed under the Boost Software License, Version 1.0.
//                 See accompanying file LICENSE.txt or copy at
//                     http://www.boost.org/LICENSE_1_0.txt
//==============================================================================
#ifndef BOOST_SIMD_REDUCTION_FUNCTIONS_SIMD_COMMON_PROD_HPP_INCLUDED
#define BOOST_SIMD_REDUCTION_FUNCTIONS_SIMD_COMMON_PROD_HPP_INCLUDED

#include <boost/simd/reduction/functions/prod.hpp>
#include <boost/simd/include/functions/simd/all_reduce.hpp>
#include <boost/simd/include/functions/simd/multiplies.hpp>

namespace boost { namespace simd { namespace ext
{
  BOOST_SIMD_FUNCTOR_IMPLEMENTATION ( boost::simd::tag::prod_, tag::cpu_
                                    , (A0)(X)
                                    , ((simd_<unspecified_<A0>,X>))
                                    )
  {
    typedef typename A0::value_type result_type;

    BOOST_FORCEINLINE result_type operator()(A0 const& a0) const
    {
      return all_reduce<boost::simd::tag::multiplies_>(a0)[0];
    }
  };
} } }

#endif

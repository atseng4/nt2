//==============================================================================
//         Copyright 2003 - 2011 LASMEA UMR 6602 CNRS/Univ. Clermont II
//         Copyright 2009 - 2011 LRI    UMR 8623 CNRS/Univ Paris Sud XI
//
//          Distributed under the Boost Software License, Version 1.0.
//                 See accompanying file LICENSE.txt or copy at
//                     http://www.boost.org/LICENSE_1_0.txt
//==============================================================================
#ifndef NT2_EXPONENTIAL_FUNCTIONS_GENERIC_LOG_HPP_INCLUDED
#define NT2_EXPONENTIAL_FUNCTIONS_GENERIC_LOG_HPP_INCLUDED

#include <nt2/exponential/functions/log.hpp>
#include <nt2/exponential/functions/scalar/impl/logs.hpp>
#include <nt2/exponential/functions/simd/common/impl/logs.hpp>
#include <nt2/include/functions/simd/tofloat.hpp>
#include <boost/simd/sdk/simd/meta/is_native.hpp>
#include <boost/dispatch/meta/as_floating.hpp>

#include <nt2/exponential/functions/generic/gruntthepeon_common.hpp>

namespace nt2 { namespace ext
{
  NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::log_, tag::cpu_
                            , (A0)
                            , (generic_< arithmetic_<A0> >)
                            )
  {
    typedef typename boost::dispatch::meta::as_floating<A0>::type result_type;
    NT2_FUNCTOR_CALL(1)
    {
      return nt2::log(nt2::tofloat(a0));
    }
  };

  NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::log_, tag::cpu_
                            , (A0)
                            , (generic_< floating_<A0> >)
                            )
  {
    typedef A0 result_type;
    typedef typename boost::simd::meta::is_native<A0>::type is_native_t;
    NT2_FUNCTOR_CALL(1)
    {
      return details::logarithm<A0,is_native_t>::log(a0);
    }
  };

  #if 1
  NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::log_, boost::simd::tag::avx_
                            , (A0)
                            , ((simd_< single_<A0>, boost::simd::tag::avx_ >))
                            )
  {
    typedef A0 result_type;
    NT2_FUNCTOR_CALL(1)
    {
      return map(functor<tag::log_>(), a0);
    }
  };

  NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::log_, boost::simd::tag::sse2_
                            , (A0)
                            , ((simd_< single_<A0>, boost::simd::tag::sse_ >))
                            )
  {
    typedef A0 result_type;
    NT2_FUNCTOR_CALL(1)
    {
      __m128 x = a0;

      __m128i emm0;
      __m128 one = *(__m128*)_ps_1;

      __m128 invalid_mask = _mm_cmple_ps(x, _mm_setzero_ps());

      x = _mm_max_ps(x, *(__m128*)_ps_min_norm_pos);  /* cut off denormalized stuff */

      emm0 = _mm_srli_epi32(_mm_castps_si128(x), 23);

      /* keep only the fractional part */
      x = _mm_and_ps(x, *(__m128*)_ps_inv_mant_mask);
      x = _mm_or_ps(x, *(__m128*)_ps_0p5);

      emm0 = _mm_sub_epi32(emm0, *(__m128i*)_pi32_0x7f);
      __m128 e = _mm_cvtepi32_ps(emm0);

      e = _mm_add_ps(e, one);

      /* part2:
         if( x < SQRTHF ) {
           e -= 1;
           x = x + x - 1.0;
         } else { x = x - 1.0; }
      */
      __m128 mask = _mm_cmplt_ps(x, *(__m128*)_ps_cephes_SQRTHF);
      __m128 tmp = _mm_and_ps(x, mask);
      x = _mm_sub_ps(x, one);
      e = _mm_sub_ps(e, _mm_and_ps(one, mask));
      x = _mm_add_ps(x, tmp);


      __m128 z = _mm_mul_ps(x,x);

      __m128 y = *(__m128*)_ps_cephes_log_p0;
      y = _mm_mul_ps(y, x);
      y = _mm_add_ps(y, *(__m128*)_ps_cephes_log_p1);
      y = _mm_mul_ps(y, x);
      y = _mm_add_ps(y, *(__m128*)_ps_cephes_log_p2);
      y = _mm_mul_ps(y, x);
      y = _mm_add_ps(y, *(__m128*)_ps_cephes_log_p3);
      y = _mm_mul_ps(y, x);
      y = _mm_add_ps(y, *(__m128*)_ps_cephes_log_p4);
      y = _mm_mul_ps(y, x);
      y = _mm_add_ps(y, *(__m128*)_ps_cephes_log_p5);
      y = _mm_mul_ps(y, x);
      y = _mm_add_ps(y, *(__m128*)_ps_cephes_log_p6);
      y = _mm_mul_ps(y, x);
      y = _mm_add_ps(y, *(__m128*)_ps_cephes_log_p7);
      y = _mm_mul_ps(y, x);
      y = _mm_add_ps(y, *(__m128*)_ps_cephes_log_p8);
      y = _mm_mul_ps(y, x);

      y = _mm_mul_ps(y, z);


      tmp = _mm_mul_ps(e, *(__m128*)_ps_cephes_log_q1);
      y = _mm_add_ps(y, tmp);


      tmp = _mm_mul_ps(z, *(__m128*)_ps_0p5);
      y = _mm_sub_ps(y, tmp);

      tmp = _mm_mul_ps(e, *(__m128*)_ps_cephes_log_q2);
      x = _mm_add_ps(x, y);
      x = _mm_add_ps(x, tmp);
      x = _mm_or_ps(x, invalid_mask); // negative arg will be NAN
      return x;
    }
  };
  #endif
} }


#endif

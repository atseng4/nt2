//==============================================================================
//         Copyright 2003 - 2011 LASMEA UMR 6602 CNRS/Univ. Clermont II
//         Copyright 2009 - 2011 LRI    UMR 8623 CNRS/Univ Paris Sud XI
//
//          Distributed under the Boost Software License, Version 1.0.
//                 See accompanying file LICENSE.txt or copy at
//                     http://www.boost.org/LICENSE_1_0.txt
//==============================================================================
#ifndef NT2_EXPONENTIAL_FUNCTIONS_GENERIC_EXP_HPP_INCLUDED
#define NT2_EXPONENTIAL_FUNCTIONS_GENERIC_EXP_HPP_INCLUDED

#include <nt2/exponential/functions/exp.hpp>
#include <nt2/exponential/functions/scalar/impl/expo.hpp>
#include <nt2/exponential/functions/simd/common/impl/expo.hpp>
#include <nt2/include/functions/simd/tofloat.hpp>
#include <boost/simd/sdk/simd/meta/is_native.hpp>
#include <boost/dispatch/meta/as_floating.hpp>

#include <nt2/exponential/functions/generic/gruntthepeon_common.hpp>

namespace nt2 { namespace ext
{
  NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::exp_, tag::cpu_
                            , (A0)
                            , (generic_< arithmetic_<A0> >)
                            )
  {
    typedef typename boost::dispatch::meta::as_floating<A0>::type result_type;
    NT2_FUNCTOR_CALL(1)
    {
      return nt2::exp(nt2::tofloat(a0));
    }
  };

  NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::exp_, tag::cpu_
                            , (A0)
                            , (generic_< floating_<A0> >)
                            )
  {
    typedef A0 result_type;
    typedef typename boost::simd::meta::is_native<A0>::type is_native_t;
    NT2_FUNCTOR_CALL(1)
    {
       return nt2::details::
              exponential<A0,natural_tag,is_native_t,accu_tag>
              ::expa(a0);
    }
  };

  #if 1
  NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::exp_, boost::simd::tag::avx_
                            , (A0)
                            , ((simd_< single_<A0>, boost::simd::tag::avx_ >))
                            )
  {
    typedef A0 result_type;
    NT2_FUNCTOR_CALL(1)
    {
      return map(functor<tag::exp_>(), a0);
    }
  };

  NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::exp_, boost::simd::tag::sse2_
                            , (A0)
                            , ((simd_< single_<A0>, boost::simd::tag::sse_ >))
                            )
  {
    typedef A0 result_type;
    NT2_FUNCTOR_CALL(1)
    {
      __m128 x = a0;

      __m128 tmp = _mm_setzero_ps(), fx;
      __m128i emm0;
      __m128 one = *(__m128*)_ps_1;

      x = _mm_min_ps(x, *(__m128*)_ps_exp_hi);
      x = _mm_max_ps(x, *(__m128*)_ps_exp_lo);

      /* express exp(x) as exp(g + n*log(2)) */
      fx = _mm_mul_ps(x, *(__m128*)_ps_cephes_LOG2EF);
      fx = _mm_add_ps(fx, *(__m128*)_ps_0p5);

      /* how to perform a floorf with SSE: just below */
      emm0 = _mm_cvttps_epi32(fx);
      tmp  = _mm_cvtepi32_ps(emm0);

      /* if greater, substract 1 */
      __m128 mask = _mm_cmpgt_ps(tmp, fx);
      mask = _mm_and_ps(mask, one);
      fx = _mm_sub_ps(tmp, mask);

      tmp = _mm_mul_ps(fx, *(__m128*)_ps_cephes_exp_C1);
      __m128 z = _mm_mul_ps(fx, *(__m128*)_ps_cephes_exp_C2);
      x = _mm_sub_ps(x, tmp);
      x = _mm_sub_ps(x, z);

      z = _mm_mul_ps(x,x);

      __m128 y = *(__m128*)_ps_cephes_exp_p0;
      y = _mm_mul_ps(y, x);
      y = _mm_add_ps(y, *(__m128*)_ps_cephes_exp_p1);
      y = _mm_mul_ps(y, x);
      y = _mm_add_ps(y, *(__m128*)_ps_cephes_exp_p2);
      y = _mm_mul_ps(y, x);
      y = _mm_add_ps(y, *(__m128*)_ps_cephes_exp_p3);
      y = _mm_mul_ps(y, x);
      y = _mm_add_ps(y, *(__m128*)_ps_cephes_exp_p4);
      y = _mm_mul_ps(y, x);
      y = _mm_add_ps(y, *(__m128*)_ps_cephes_exp_p5);
      y = _mm_mul_ps(y, z);
      y = _mm_add_ps(y, x);
      y = _mm_add_ps(y, one);

      /* build 2^n */
      emm0 = _mm_cvttps_epi32(fx);
      emm0 = _mm_add_epi32(emm0, *(__m128i*)_pi32_0x7f);
      emm0 = _mm_slli_epi32(emm0, 23);
      __m128 pow2n = _mm_castsi128_ps(emm0);
      y = _mm_mul_ps(y, pow2n);
      return y;
    }
  };
  #endif
} }


#endif

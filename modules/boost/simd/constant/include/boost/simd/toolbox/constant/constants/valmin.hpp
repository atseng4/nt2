//==============================================================================
//         Copyright 2003 - 2011   LASMEA UMR 6602 CNRS/Univ. Clermont II
//         Copyright 2009 - 2011   LRI    UMR 8623 CNRS/Univ Paris Sud XI
//
//          Distributed under the Boost Software License, Version 1.0.
//                 See accompanying file LICENSE.txt or copy at
//                     http://www.boost.org/LICENSE_1_0.txt
//==============================================================================
/*!
 * \file
**/
#ifndef BOOST_SIMD_TOOLBOX_CONSTANT_CONSTANTS_VALMIN_HPP_INCLUDED
#define BOOST_SIMD_TOOLBOX_CONSTANT_CONSTANTS_VALMIN_HPP_INCLUDED

#include <boost/simd/include/simd.hpp>
#include <boost/simd/sdk/meta/int_c.hpp>
#include <boost/simd/sdk/meta/float.hpp>
#include <boost/simd/sdk/meta/double.hpp>
#include <boost/simd/sdk/constant/common.hpp>
#include <boost/simd/sdk/constant/constant.hpp>

#ifdef BOOST_MSVC
#pragma warning(push)
#pragma warning(disable: 4146)
#endif

/*!
 * \ingroup boost_simd_constant
 * \defgroup boost_simd_constant_valmin Valmin
 *
 * \par Description
 * Constant Valmin, maximum value of a type.
 * \arg int8    -128, uint8    0,
 * \arg int16 -32768, uint16 0,
 * \arg int32 -2147483648, uint32 0,
 * \arg int64 -9223372036854775808, uint64 0,\arg float \f$-\infty\f$, double \f$-\infty\f$,
 * \par
 * The value of this constant is type dependant. This means that for different
 * types it does not represent the same mathematical number.
 *
 * \par Header file
 * 
 * \code
 * #include <nt2/include/functions/valmin.hpp>
 * \endcode
 * 
 * 
 * \synopsis
 *
 * \code
 * namespace boost::simd
 * {
 *   template <class T,class A0>
 *     meta::call<tag::valmin_(A0)>::type
 *     Valmin();
 * }
 * \endcode
 *
 * 
 * \param T template parameter of Valmin
 * 
 * \return type T value
 *  
 *  
**/

namespace boost { namespace simd
{
  namespace tag
  {
    /*!
     * \brief Define the tag Valmin of functor Valmin 
     *        in namespace boost::simd::tag for toolbox boost.simd.constant
    **/
    struct Valmin
    { 
      typedef double default_type;
      template<class Target, class Dummy=void> 
      struct apply : meta::int_c<Target,0> {}; 
    };

    template<class Dummy>
    struct  Valmin::apply<float,Dummy> 
          : meta::single_<0xFF7FFFFFUL> {};

    template<class Dummy>
    struct  Valmin::apply<double,Dummy> 
          : meta::double_<0xFFEFFFFFFFFFFFFFULL> {};

    template<class Dummy>
    struct  Valmin::apply<boost::simd::int8_t,Dummy> 
          : meta::int_c<boost::simd::int8_t,-boost::simd::uint8_t(128)> {};

    template<class Dummy>
    struct  Valmin::apply<boost::simd::int16_t,Dummy> 
          : meta::int_c<boost::simd::int16_t,-boost::simd::uint16_t(32768U)> {};

    template<class Dummy>
    struct  Valmin::apply<boost::simd::int32_t,Dummy> 
          : meta::int_c<boost::simd::int32_t,-boost::simd::uint32_t(2147483648UL)> {};

    template<class Dummy>
    struct  Valmin::apply<boost::simd::int64_t,Dummy> 
          : meta::int_c<boost::simd::int64_t,-boost::simd::uint64_t(9223372036854775808ULL)> {};
  }

  BOOST_SIMD_CONSTANT_IMPLEMENTATION(boost::simd::tag::Valmin, Valmin)
} }

#ifdef BOOST_MSVC
#pragma warning(pop)
#endif

#endif

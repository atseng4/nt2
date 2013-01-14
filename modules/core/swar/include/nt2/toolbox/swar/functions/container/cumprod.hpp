//==============================================================================
//         Copyright 2003 - 2012   LASMEA UMR 6602 CNRS/Univ. Clermont II
//         Copyright 2009 - 2012   LRI    UMR 8623 CNRS/Univ Paris Sud XI
//
//          Distributed under the Boost Software License, Version 1.0.
//                 See accompanying file LICENSE.txt or copy at
//                     http://www.boost.org/LICENSE_1_0.txt
//==============================================================================
#ifndef NT2_TOOLBOX_SWAR_FUNCTIONS_CONTAINER_CUMPROD_HPP_INCLUDED
#define NT2_TOOLBOX_SWAR_FUNCTIONS_CONTAINER_CUMPROD_HPP_INCLUDED

#include <nt2/toolbox/swar/functions/cumprod.hpp>
#include <boost/simd/toolbox/swar/functions/cumprod.hpp>
#include <nt2/sdk/meta/size_as.hpp>
#include <nt2/sdk/meta/value_as.hpp>
#include <nt2/core/container/dsl/size.hpp>
#include <nt2/core/container/dsl/value_type.hpp>

namespace nt2 { namespace ext
{
  template<class Domain, int N, class Expr>
  struct size_of<boost::simd::tag::cumprod_,Domain, N,Expr>
        : meta::size_as<Expr,0>
  {};

  template<class Domain, int N, class Expr>
  struct  value_type<boost::simd::tag::cumprod_,Domain,N,Expr>
        : meta::value_as<Expr,0>
  {};
} }

#endif

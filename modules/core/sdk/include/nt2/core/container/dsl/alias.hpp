//==============================================================================
//         Copyright 2003 - 2011   LASMEA UMR 6602 CNRS/Univ. Clermont II
//         Copyright 2009 - 2011   LRI    UMR 8623 CNRS/Univ Paris Sud XI
//         Copyright 2012 - 2013   MetaScale
//
//          Distributed under the Boost Software License, Version 1.0.
//                 See accompanying file LICENSE.txt or copy at
//                     http://www.boost.org/LICENSE_1_0.txt
//==============================================================================
#ifndef NT2_CORE_CONTAINER_DSL_ALIAS_HPP_INCLUDED
#define NT2_CORE_CONTAINER_DSL_ALIAS_HPP_INCLUDED

#include <nt2/include/functions/numel.hpp>
#include <nt2/core/container/dsl/expression.hpp>

namespace nt2 { namespace container
{
  template<class T, class U>
  bool alias(T const& t, U const& u)
  {
    return raw(value(t)) < raw(value(u))+numel(u)
        && raw(value(t))+numel(t) >= raw(value(u))
    ;
  }
} }

#endif

/*******************************************************************************
 *         Copyright 2003-2010 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2009-2010 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_SDK_DSL_CALL_HPP_INCLUDED
#define NT2_SDK_DSL_CALL_HPP_INCLUDED

////////////////////////////////////////////////////////////////////////////////
// This file generate basic EDSL expression wrapper over any nt2 function
////////////////////////////////////////////////////////////////////////////////
#include <nt2/sdk/meta/any.hpp>
#include <boost/proto/proto.hpp>
#include <nt2/sdk/dsl/category.hpp>
#include <nt2/sdk/functor/functor.hpp>
#include <nt2/sdk/functor/meta/call.hpp>
#include <nt2/sdk/dsl/proto/as_child.hpp>
#include <nt2/sdk/functor/meta/hierarchy.hpp>

#if defined(NT2_DONT_USE_PREPROCESSED_FILES)
#include <boost/dispatch/extension/parameters.hpp>
#include <boost/preprocessor/selection/min.hpp>
#include <nt2/sdk/functor/preprocessor/call.hpp>
#include <nt2/sdk/functor/preprocessor/dispatch.hpp>
#include <boost/preprocessor/repetition/repeat_from_to.hpp>
#include <boost/preprocessor/repetition/enum_binary_params.hpp>
#endif

////////////////////////////////////////////////////////////////////////////////
// Defines the catch-all call for proto expression
////////////////////////////////////////////////////////////////////////////////
#if !defined(NT2_DONT_USE_PREPROCESSED_FILES)
#include <nt2/sdk/dsl/preprocessed/call.hpp>
#else
#if defined(__WAVE__) && defined(NT2_CREATE_PREPROCESSED_FILES) && __INCLUDE_LEVEL__ == 0
#pragma wave option(preserve: 2, line: 0, output: "preprocessed/call.hpp")
#endif

#define M1(z,n,t) nt2::meta::as_child(BOOST_PP_CAT(a,n))
#define M2(z,n,t) (BOOST_PP_CAT(A,n))
#define M3(z,n,t) (unspecified_<BOOST_PP_CAT(A,n)>)

#define M4(z,n,t)                                                             \
NT2_REGISTER_DISPATCH_IF( Func, tag::formal_                                  \
                        , (Func)BOOST_PP_REPEAT(n,M2,~)                       \
                        , (any< boost::proto::is_expr<boost::mpl::_>          \
                              , BOOST_PP_ENUM_PARAMS(n,A)                     \
                             >                                                \
                          )                                                   \
                      , (Func(tag::ast_))                                     \
                      , BOOST_PP_REPEAT(n,M3,~)                               \
                      )                                                       \
/**/

#define M0(z,n,t)                                                 \
template<class This,BOOST_PP_ENUM_PARAMS(n,class A)>              \
struct result<This(BOOST_PP_ENUM_PARAMS(n,A))>                    \
{                                                                 \
  typedef typename boost::proto::result_of::                      \
  make_expr < Func                                                \
            , BOOST_PP_ENUM_BINARY_PARAMS                         \
              ( n                                                 \
              , typename nt2::details::result_of                  \
                ::as_child< typename meta::strip< A               \
              ,                                 >::type const&    \
                          >::type BOOST_PP_INTERCEPT              \
              )                                                   \
            >::type type;                                         \
};                                                                \
template<BOOST_PP_ENUM_PARAMS(n,class A)> inline                  \
typename result<implement                                         \
                (BOOST_PP_ENUM_BINARY_PARAMS( n,A                 \
                                            , const&              \
                                              BOOST_PP_INTERCEPT) \
                )                                                 \
               >::type                                            \
operator()(BOOST_PP_ENUM_BINARY_PARAMS(n,A,const& a) ) const      \
{                                                                 \
  return boost::proto::                                           \
  make_expr<Func>( BOOST_PP_ENUM(n,M1,~) );                       \
}                                                                 \
/**/

namespace nt2 { namespace meta
{
  BOOST_PP_REPEAT_FROM_TO ( 1
                          , BOOST_PP_INC(BOOST_PP_MIN ( NT2_MAX_ARITY
                                                      , BOOST_PROTO_MAX_ARITY
                                                      )
                                        )
                         ,M4,~
                         )

  template<class Func,class Dummy>
  struct implement<Func(tag::ast_),tag::formal_,Dummy>
  {
    template<class Sig> struct result;
    BOOST_PP_REPEAT_FROM_TO ( 1
                            , BOOST_PP_INC(BOOST_PP_MIN ( NT2_MAX_ARITY
                                                        , BOOST_PROTO_MAX_ARITY
                                                        )
                                          )
                           ,M0,~
                           )
  };
} }

#undef M0
#undef M1
#undef M2
#undef M3
#undef M4

#if defined(__WAVE__) && defined(NT2_CREATE_PREPROCESSED_FILES)
#pragma wave option(output: null)
#endif
#endif

#endif


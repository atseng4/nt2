//==============================================================================
//         Copyright 2003 - 2012   LASMEA UMR 6602 CNRS/Univ. Clermont II
//         Copyright 2009 - 2013   LRI    UMR 8623 CNRS/Univ Paris Sud XI
//
//          Distributed under the Boost Software License, Version 1.0.
//                 See accompanying file LICENSE.txt or copy at
//                     http://www.boost.org/LICENSE_1_0.txt
//==============================================================================
#ifndef NT2_TOOLBOX_LINALG_FUNCTIONS_LAPACK_POSV_HPP_INCLUDED
#define NT2_TOOLBOX_LINALG_FUNCTIONS_LAPACK_POSV_HPP_INCLUDED

#include <nt2/linalg/functions/posv.hpp>
#include <nt2/include/functions/zeros.hpp>
#include <nt2/core/container/table/category.hpp>
#include <nt2/linalg/details/utility/workspace.hpp>
#include <nt2/dsl/functions/terminal.hpp>
#include <nt2/include/functions/width.hpp>
#include <nt2/linalg/details/utility/f77_wrapper.hpp>

#include <nt2/table.hpp>

extern "C"
{
  void NT2_F77NAME(dposv)(  const char* uplo         , const nt2_la_int* n
                          , const nt2_la_int* nhrs   , double* a
                          , const nt2_la_int* lda    , const double* b
                          , const nt2_la_int* ldb    , nt2_la_int* info
                          );

  void NT2_F77NAME(sposv)(  const char* uplo         , const nt2_la_int* n
                          , const nt2_la_int* nhrs   , float* a
                          , const nt2_la_int* lda    , const float* b
                          , const nt2_la_int* ldb    , nt2_la_int* info
                          );
}

namespace nt2 { namespace ext
{
  /// INTERNAL ONLY - Compute the workspace
  NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::posv_, tag::cpu_
                            , (A0)(S0)(A1)(S1)(A2)(S2)
                            , ((expr_ < table_< double_<A0>, S0 >     // A
                                      , nt2::tag::terminal_
                                      , boost::mpl::long_<0>
                                      >
                              ))
                              ((expr_ < table_< double_<A1>, S1 >     // B
                                      , nt2::tag::terminal_
                                      , boost::mpl::long_<0>
                                      >
                              ))
                              ((expr_ < table_< double_<A2>, S2 >     // X
                                      , nt2::tag::terminal_
                                      , boost::mpl::long_<0>
                                      >
                              ))
                            )
  {
     typedef nt2_la_int result_type;

     BOOST_FORCEINLINE result_type operator()(A0& a0, A1 const& a1, A2& a2) const
     {
        result_type that;
        nt2_la_int  n  = nt2::width(a0);
        nt2_la_int  ld = a0.leading_size();
        nt2_la_int  ldb = a1.leading_size();
        nt2_la_int  nhrs = nt2::width(a1);
        char uplo = 'L';

        NT2_F77NAME(dposv) ( &uplo, &n, &nhrs, a0.raw(), &ld, a2.raw(), &ldb
                           , &that
                            );

        return that;
     }
  };



  /// INTERNAL ONLY - Compute the workspace
  NT2_FUNCTOR_IMPLEMENTATION( nt2::tag::posv_, tag::cpu_
                            , (A0)(S0)(A1)(S1)(A2)(S2)
                            , ((expr_ < table_< single_<A0>, S0 >     // A
                                      , nt2::tag::terminal_
                                      , boost::mpl::long_<0>
                                      >
                              ))
                              ((expr_ < table_< single_<A1>, S1 >     // B
                                      , nt2::tag::terminal_
                                      , boost::mpl::long_<0>
                                      >
                              ))
                              ((expr_ < table_< single_<A2>, S2 >     // X
                                      , nt2::tag::terminal_
                                      , boost::mpl::long_<0>
                                      >
                              ))
                            )
  {
     typedef nt2_la_int result_type;

     BOOST_FORCEINLINE result_type operator()(A0& a0, A1 const& a1, A2& a2) const
     {
        result_type that;
        details::workspace<typename A0::value_type> w;
        nt2_la_int  n  = nt2::width(a0);
        nt2_la_int  ld = a0.leading_size();
        nt2_la_int  ldb = a1.leading_size();
        nt2_la_int  nhrs = nt2::width(a1);
        char uplo = 'L';


        NT2_F77NAME(sposv) ( &uplo, &n, &nhrs, a0.raw(), &ld, a2.raw(), &ldb
                           , &that
                            );

        return that;
     }
  };


} }

#endif

//==============================================================================
//         Copyright 2003 - 2011   LASMEA UMR 6602 CNRS/Univ. Clermont II
//         Copyright 2009 - 2013   LRI    UMR 8623 CNRS/Univ Paris Sud XI
//         Copyright 2012 - 2013   MetaScale SAS
//
//          Distributed under the Boost Software License, Version 1.0.
//                 See accompanying file LICENSE.txt or copy at
//                     http://www.boost.org/LICENSE_1_0.txt
//==============================================================================
#if !BOOST_PP_IS_ITERATING
#ifndef NT2_SDK_TBB_FUTURE_WHEN_ALL_HPP_INCLUDED
#define NT2_SDK_TBB_FUTURE_WHEN_ALL_HPP_INCLUDED

#if defined(NT2_USE_TBB)

#include <tbb/tbb.h>
#include <tbb/flow_graph.h>

#include <vector>
#include <nt2/sdk/tbb/future/details/tbb_future.hpp>

namespace nt2
{
  namespace details
  {

    struct empty_body
    {
      int operator()(){}
    };

    template<class Site>
    struct when_all_impl< tag::tbb_<Site> >
    {
       typedef typename tbb::flow::continue_node<\
       tbb::flow::continue_msg> node_type;

#define BOOST_PP_ITERATION_PARAMS_1 (3, \
( 1, BOOST_DISPATCH_MAX_ARITY, \
"nt2/sdk/tbb/future/details/tbb_when_all.hpp") \
)

#include BOOST_PP_ITERATE()
    };
  }
}

#endif
#endif

#else

#define N BOOST_PP_ITERATION()

#define NT2_FUTURE_FORWARD_ARGS(z,n,t) details::tbb_future<A##n> const & a##n
#define NT2_FUTURE_FORWARD_ARGS1(z,n,t) tbb::flow::make_edge(*a##n##.node_,*c);

        template< BOOST_PP_ENUM_PARAMS(N, typename A) >
        details::tbb_future<int> call\
        ( BOOST_PP_ENUM(N,NT2_FUTURE_FORWARD_ARGS, ~))
        {
            tbb_future<int> future_res;

            details::empty_body f;

            tbb::flow::graph * work = a0.work_;

            std::vector<node_type *> * node_list \
              = a0.node_list_;

            bool * ready = a0.ready_;

            node_type * c = new node_type
                ( *work_,
                  details::tbb_task_wrapper0<F,int>
                  (f,future_res.res_)
                );

            node_list_->push_back(c);

            BOOST_PP_REPEAT(N, NT2_FUTURE_FORWARD_ARGS2, ~)

            future_res.attach_task(work,node_list,c,ready);

            return future_res;
         }

#undef NT2_FUTURE_FORWARD_ARGS
#undef NT2_FUTURE_FORWARD_ARGS1

#endif

//==============================================================================
//         Copyright 2003 - 2011   LASMEA UMR 6602 CNRS/Univ. Clermont II
//         Copyright 2009 - 2013   LRI    UMR 8623 CNRS/Univ Paris Sud XI
//         Copyright 2012 - 2013   MetaScale SAS
//
//          Distributed under the Boost Software License, Version 1.0.
//                 See accompanying file LICENSE.txt or copy at
//                     http://www.boost.org/LICENSE_1_0.txt
//==============================================================================
#ifndef NT2_SDK_TBB_SPAWNER_TRANSFORM_HPP_INCLUDED
#define NT2_SDK_TBB_SPAWNER_TRANSFORM_HPP_INCLUDED

#if defined(NT2_USE_TBB)

#include <tbb/tbb.h>
#include <nt2/sdk/tbb/blocked_range.hpp>
#include <nt2/sdk/shared_memory/spawner.hpp>

#include <nt2/sdk/tbb/future/future.hpp>
#include <cstdio>

#ifndef BOOST_NO_EXCEPTIONS
#include <boost/exception_ptr.hpp>
#endif

struct p1
{
    typedef int result_type;

    int operator()()
    {
        return 1;
    }
};

struct p2
{
    typedef int result_type;

    int operator()(int value)
    {
        return value*2;
    }
};

struct p3
{
    typedef int result_type;

    int operator()(int value)
    {
        printf("Result is: %d\n",value);
        return value;
    }
};


namespace nt2
{
    namespace tag
    {
        struct transform_;
        template<class T> struct tbb_;
    }

//    namespace details
//    {
//        template<class Worker>
//        struct Tbb_Transformer
//        {
//            Tbb_Transformer(Worker & w)
//            :w_(w)
//            {}
//
//            void operator()(nt2::blocked_range<std::size_t> const& r) const
//            {
//                w_(r.begin(),r.size());
//            };
//
//            Worker & w_;
//
//            private:
//            Tbb_Transformer& operator=(Tbb_Transformer const&);
//
//        };
//    }

    template<class Site>
    struct spawner< tag::transform_, tag::tbb_<Site> >
    {
        typedef typename tag::tbb_<Site> Arch;

        spawner(){}

        template<typename Worker>
        void operator()(Worker & w, std::size_t begin, std::size_t size, std::size_t grain)
        {

//#ifndef BOOST_NO_EXCEPTIONS
//            boost::exception_ptr exception;
//
//            try
//            {
//#endif
//            details::Tbb_Transformer<Worker> tbb_w ( w );
//
//            tbb::parallel_for(
//                nt2::blocked_range<std::size_t>(begin,begin+size,grain),
//                tbb_w);
//
//#ifndef BOOST_NO_EXCEPTIONS
//            }
//            catch(...)
//            {
//            exception = boost::current_exception();
//            }
//#endif

        typedef typename
        nt2::details::tbb_future<int> future;

        p1 w1;
        p2 w2;

        future f1 = async<Arch>(w1);
        future f2 = f1.then(w2);

        f2.get();

        }
    };
}

#endif
#endif

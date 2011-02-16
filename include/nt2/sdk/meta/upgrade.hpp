/*******************************************************************************
 *         Copyright 2003 & onward LASMEA UMR 6602 CNRS/Univ. Clermont II
 *         Copyright 2009 & onward LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
#ifndef NT2_SDK_META_UPGRADE_HPP_INCLUDED
#define NT2_SDK_META_UPGRADE_HPP_INCLUDED

//////////////////////////////////////////////////////////////////////////////
// For any type, return the integer type of size equals to sizeof(T)*2
// with an optional sign change, or the real type if input is real
// See: http://nt2.metascale.org/sdk/meta/traits/upgrade.html
//////////////////////////////////////////////////////////////////////////////
#include <nt2/sdk/meta/strip.hpp>
#include <nt2/sdk/meta/sign_of.hpp>
#include <nt2/sdk/meta/make_integer.hpp>
#include <nt2/sdk/meta/primitive_of.hpp>
#include <nt2/sdk/meta/factory_of.hpp>
#include <nt2/sdk/meta/is_unspecified.hpp>
#include <boost/mpl/apply.hpp>

namespace nt2 { namespace details
{
  //////////////////////////////////////////////////////////////////////////////
  // Integral types are upgraded using make_integer
  //////////////////////////////////////////////////////////////////////////////
  template<class T,std::size_t Size, class Sign, class Lambda>
  struct upgrade : meta::make_integer<Size*2,Sign,Lambda> {};

  //////////////////////////////////////////////////////////////////////////////
  // If type size is 8, return the type itself for any category
  //////////////////////////////////////////////////////////////////////////////
  template<class T, class Sign, class Lambda>
  struct   upgrade<T,8,Sign,Lambda>
        : meta::make_integer<8,Sign,Lambda> {};


  template<class Sign, class Lambda>
  struct upgrade<double,sizeof(double),Sign,Lambda>
  {
    typedef typename boost::mpl::apply1<Lambda, double>::type type;
  };

  template<class Sign, class Lambda>
  struct upgrade<float,sizeof(float),Sign,Lambda>
  {
    typedef typename boost::mpl::apply1<Lambda, double>::type type;
  };
} }

namespace nt2 { namespace meta
{
  //////////////////////////////////////////////////////////////////////////////
  // For any type, return the integer type of size equals to sizeof(T)*2
  // with an optional sign change
  //////////////////////////////////////////////////////////////////////////////
  template<class T,class Sign=typename meta::sign_of<T>::type>
  struct  upgrade
        : details::upgrade< typename meta::primitive_of<typename meta::strip<T>::type>::type
                          , sizeof(typename meta::primitive_of<typename meta::strip<T>::type>::type)
                          ,Sign
                          , typename meta::factory_of<typename meta::strip<T>::type>::type
                          >
  {
    NT2_STATIC_ASSERT ( (!is_unspecified<T>::value)
                      , NT2_UNHIERARCHIZED_TYPE_USED_IN_META_UPGRADE
                      , "An unhierarchized type is used in nt2::meta::upgrade."
                      );
  };

} }

#endif


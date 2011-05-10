//////////////////////////////////////////////////////////////////////////////
///   Copyright 2003 and onward LASMEA UMR 6602 CNRS/U.B.P Clermont-Ferrand
///   Copyright 2009 and onward LRI    UMR 8623 CNRS/Univ Paris Sud XI
///
///          Distributed under the Boost Software License, Version 1.0
///                 See accompanying file LICENSE.txt or copy at
///                     http://www.boost.org/LICENSE_1_0.txt
//////////////////////////////////////////////////////////////////////////////
#ifndef NT2_TOOLBOX_DUAL_FUNCTION_SCALAR_COMPLEMENT_HPP_INCLUDED
#define NT2_TOOLBOX_DUAL_FUNCTION_SCALAR_COMPLEMENT_HPP_INCLUDED
/////////////////////////////////////////////////////////////////////////////
// Implementation when type A0 is dual
/////////////////////////////////////////////////////////////////////////////
NT2_REGISTER_DISPATCH(tag::complement_, tag::cpu_,
                      (A0),
                      ((dual_<real_<A0> > ))
                     ); 

namespace nt2 { namespace ext
{
  template<class Dummy>
  struct call<tag::complement_(tag::dual_ <tag::real_>),
              tag::cpu_, Dummy> : callable
	      {
    template<class Sig> struct result;
    template<class This,class A0>
    struct result<This(A0)> : meta::strip<A0>{};

    NT2_FUNCTOR_CALL(1)
    {
      A0 t(complement(nt2::get<0>(a0)), complement(nt2::get<1>(a0))); 
      return t; 
    }
  };
} }

#endif

//////////////////////////////////////////////////////////////////////////////
///   Copyright 2003 or onward LASMEA UMR 6602 CNRS/U.B.P Clermont-Ferror
///   Copyright 2009 or onward LRI    UMR 8623 CNRS/Univ Paris Sud XI
///
///          Distributed under the Boost Software License, Version 1.0
///                 See accompanying file LICENSE.txt or copy at
///                     http://www.boost.org/LICENSE_1_0.txt
//////////////////////////////////////////////////////////////////////////////
#ifndef NT2_TOOLBOX_DUAL_FUNCTION_SCALAR_BITWISE_OR_HPP_INCLUDED
#define NT2_TOOLBOX_DUAL_FUNCTION_SCALAR_BITWISE_OR_HPP_INCLUDED

/////////////////////////////////////////////////////////////////////////////
// Implementation when type A0 or A1 are dual
/////////////////////////////////////////////////////////////////////////////
NT2_REGISTER_DISPATCH(tag::bitwise_or_, tag::cpu_,
                      (A0),
                      ((dual_<real_<A0> > ))
                      ((dual_<real_<A0> > )) 
                     ); 

namespace nt2 { namespace ext
{
  template<class Dummy>
  struct call<tag::bitwise_or_(tag::dual_ <tag::real_> ,
			tag::dual_ <tag::real_> ),
              tag::cpu_, Dummy> : callable
  {
    template<class Sig> struct result;
    template<class This,class A0>
    struct result<This(A0, A0)> : meta::strip<A0>{};

    NT2_FUNCTOR_CALL(2)
    {
      return A0(b_or(nt2::get<0>(a0), nt2::get<0>(a1)),
		b_or(nt2::get<1>(a0), nt2::get<1>(a1))); 
    }
  };
} }
/////////////////////////////////////////////////////////////////////////////
// Implementation when type A0 is dual
/////////////////////////////////////////////////////////////////////////////
NT2_REGISTER_DISPATCH(tag::bitwise_or_, tag::cpu_,
                      (A0),
                      ((dual_<real_<A0> > ))
                      ((real_<A0> )) 
                     ); 

namespace nt2 { namespace ext
{
  template<class Dummy>
  struct call<tag::bitwise_or_(tag::dual_ <tag::real_> ,
			 tag::real_ ),
              tag::cpu_, Dummy> : callable
  {
    template<class Sig> struct result;
    template<class This,class A0, class A1>
    struct result<This(A0, A1)> : meta::strip<A0>{};

    NT2_FUNCTOR_CALL(2)
    {
      return A0(b_or(nt2::get<0>(a0), a1),
		nt2::get<1>(a0)); 
    }
  };
} }
/////////////////////////////////////////////////////////////////////////////
// Implementation when type A1 is dual
/////////////////////////////////////////////////////////////////////////////
NT2_REGISTER_DISPATCH(tag::bitwise_or_, tag::cpu_,
                      (A0),
                      ((real_<A0> )) 
                      ((dual_<real_<A0> > ))
                     ); 

namespace nt2 { namespace ext
{
  template<class Dummy>
  struct call<tag::bitwise_or_(tag::real_,
			 tag::dual_ <tag::real_> 
			 ),
              tag::cpu_, Dummy> : callable
  {
    template<class Sig> struct result;
    template<class This,class A0, class A1>
    struct result<This(A0, A1)> : meta::strip<A1>{};

    NT2_FUNCTOR_CALL(2)
    {
      return A0(b_or(nt2::get<0>(a1), a0),
		nt2::get<1>(a1)); 
    }
  };
} }


#endif

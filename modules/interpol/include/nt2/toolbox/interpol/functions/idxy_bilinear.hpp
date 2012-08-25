/*******************************************************************************
 *         Copyright 2003-2012 LASMEA UMR 6602 CNRS/U.B.P
 *         Copyright 2011-2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
 *
 *          Distributed under the Boost Software License, Version 1.0.
 *                 See accompanying file LICENSE.txt or copy at
 *                     http://www.boost.org/LICENSE_1_0.txt
 ******************************************************************************/
/*!
 * \file
**/
#ifndef NT2_TOOLBOX_STATISTICS_FUNCTIONS_IDXY_BILINEAR_HPP_INCLUDED
#define NT2_TOOLBOX_STATISTICS_FUNCTIONS_IDXY_BILINEAR_HPP_INCLUDED
#include <nt2/include/simd.hpp>
#include <nt2/include/functor.hpp>
#include <boost/mpl/bool.hpp>
#include <nt2/sdk/meta/adapted_traits.hpp>

/*!
 * \ingroup interpol
 * \defgroupinterpol_idxy_bilinear idxy_bilinear
 *
 * \par Description
 * one dimensional idxy_bilinear interpolation
 * given an array a and a vector idx of "real" indices
 * provides the interpolated values along the a rows
 * by bilinear formula.
 * calls can be idxy_bilinear(a0, idx, idy)
 *              idxy_bilinear(a0, idx, idy, true) allowing extrapolation
 *              idxy_bilinear(a0, idx, idy, val1) putting val1 outside the bounds
 *              idxy_bilinear(a0, idx, idy, val1, val2)}, putting val1 under the index bounds val2 over
 *              idxy_bilinear(a0, idx, idy, val1x val2x, val1y, val2y)
 *              idxy_bilinear(a0, idx, idy, _, dim1, dim2) are the dimensions of interpolation 2 (rows), 1 (columns) of a0 by default
 *              idxy_bilinear(a0, idx, idy, val1, dim1, dim2)
 *              idxy_bilinear(a0, idx, idy, val1, val2, dim1, dim2)
 *              idxy_bilinear(a0, idx, idy, val1x val2x, val1y, val2y, dim1, dim2) can also be used
 * \par
 *
 * \par Header file
 * 
 * \code
 * #include <nt2/include/functions/idxy_bilinear.hpp>
 * \endcode
 * 
 *  
**/

namespace nt2 { namespace tag
  {         
    /*!
     * \brief Define the tag idxy_bilinear_ of functor idxy_bilinear 
     *        in namespace nt2::tag for toolbox statistics
    **/
    struct idxy_bilinear_ : ext::unspecified_<idxy_bilinear_> { typedef ext::unspecified_<idxy_bilinear_> parent; };
  }
  NT2_FUNCTION_IMPLEMENTATION(tag::idxy_bilinear_, idxy_bilinear, 3)  
  NT2_FUNCTION_IMPLEMENTATION(tag::idxy_bilinear_, idxy_bilinear, 4)
  NT2_FUNCTION_IMPLEMENTATION(tag::idxy_bilinear_, idxy_bilinear, 5)
  NT2_FUNCTION_IMPLEMENTATION(tag::idxy_bilinear_, idxy_bilinear, 6)  
  NT2_FUNCTION_IMPLEMENTATION(tag::idxy_bilinear_, idxy_bilinear, 7)
  NT2_FUNCTION_IMPLEMENTATION(tag::idxy_bilinear_, idxy_bilinear, 9)
  NT2_FUNCTION_IMPLEMENTATION(tag::idxy_bilinear_, idxy_bilinear, 10)
}

namespace nt2 { namespace container { namespace ext
{
  template<class Domain, class Expr,  int N>
  struct size_of<tag::idxy_bilinear_, Domain, N, Expr>
  {
    typedef typename boost::proto::result_of::child_c<Expr&,0>::type         value_t;
    typedef typename boost::proto::result_of::child_c<Expr&,1>::type           idx_t;
    typedef typename boost::proto::result_of::child_c<Expr&,2>::type           idy_t;
    typedef typename boost::proto::result_of::child_c<Expr&,N-2>::type   choice1_t; 
    typedef typename boost::proto::result_of::child_c<Expr&,N-1>::type   choice2_t; 

    typedef typename boost::proto::result_of::child_c<Expr&, 1>::value_type child1;
    typedef typename child1::extent_type                               result_type; 

    result_type operator()(Expr& e) const
    {
      const value_t & y     =  boost::proto::child_c<0>(e);
      const idx_t & xi      =  boost::proto::child_c<1>(e);
      const idy_t & yi      =  boost::proto::child_c<2>(e);
      const choice1_t & v1   =  boost::proto::child_c<N-2>(e);
      const choice2_t & v2   =  boost::proto::child_c<N-1>(e);
      typedef typename nt2::meta::is_integral<choice1_t>::type c_t; 
      std::size_t dim1 = 2, dim2 = 1;
      getdims(xi, v1, v2, dim1, dim2, c_t());
      result_type sizee = y.extent(); 
      sizee[dim1-1] = numel(xi);
      sizee[dim2-1] = numel(yi);
      return sizee;
    }
    static void getdims(const idx_t & xi,
                        const choice1_t &, const choice2_t &,
                        size_t&, size_t&,
                        const boost::mpl::false_& )
    { }
    static void getdims(const idx_t &,
                        const choice1_t &v1, const choice2_t &v2,
                        size_t& dim1, size_t& dim2, 
                        const boost::mpl::true_& )
    {dim1 = v1; dim2 = v2; }
  
    
  };

 template <class Domain, class Expr,  int N>
 struct value_type < tag::idxy_bilinear_, Domain,N,Expr> {
   typedef typename boost::proto::result_of::child_c<Expr&, 0>::value_type  child1;
   typedef typename nt2::meta::scalar_of<child1>::type                    elt_type;
   typedef typename nt2::meta::strip<elt_type>::type                          type;
 }; 
} } }
#endif

// /////////////////////////////////////////////////////////////////////////////
// End of idxy_bilinear.hpp
// /////////////////////////////////////////////////////////////////////////////
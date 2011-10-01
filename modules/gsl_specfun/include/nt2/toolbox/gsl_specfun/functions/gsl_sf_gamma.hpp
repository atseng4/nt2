//////////////////////////////////////////////////////////////////////////////
///   Copyright 2003 and onward LASMEA UMR 6602 CNRS/U.B.P Clermont-Ferrand
///   Copyright 2009 and onward LRI    UMR 8623 CNRS/Univ Paris Sud XI
///
///          Distributed under the Boost Software License, Version 1.0
///                 See accompanying file LICENSE.txt or copy at
///                     http://www.boost.org/LICENSE_1_0.txt
//////////////////////////////////////////////////////////////////////////////
/*!
 * \file
**/
#ifndef NT2_TOOLBOX_GSL_SPECFUN_FUNCTIONS_GSL_SF_GAMMA_HPP_INCLUDED
#define NT2_TOOLBOX_GSL_SPECFUN_FUNCTIONS_GSL_SF_GAMMA_HPP_INCLUDED
#include <nt2/include/simd.hpp>
#include <nt2/include/functor.hpp>

/*!
 * \ingroup gsl_specfun
 * \defgroup gsl_specfun_gsl_sf_gamma gsl_sf_gamma function
 *
 * \par Description
 * TODO Put description here
 *
 * \par Header file
 * 
 * \code
 * #include <nt2/toolbox/gsl_specfun/include/functions/gsl_sf_gamma.hpp>
 * \endcode
 * 
 * 
 * \synopsis
 *
 * \code
 * namespace nt2
 * {
 *   namespace gsl_specfun
 *   {
 *     template <class A0>
 *       meta::call<tag::gsl_sf_gamma_(A0)>::type
 *       gsl_sf_gamma(const A0 & a0);
 *   }
 * }
 * \endcode
 *
 * \param a0 the unique parameter of gsl_sf_gamma
 * 
 * \return a value of the same type as the parameter
 *  
 * \par Notes
 * In SIMD mode, this function acts elementwise on the inputs vectors elements
 * \par
 * \par
 * gsl_specfun library defines functions for float and double entries.
 * \par
 * As they are written in C the original name of the float version is
 * generally terminated by and extra 'f',
 * this is not the case for the nt2 version which dispatch to
 * the correct function according to the inputs types.
 *  
**/

namespace nt2 { namespace gsl_specfun { namespace tag
  {         
    /*!
     * \brief Define the tag gsl_sf_gamma_ of functor gsl_sf_gamma 
     *        in namespace nt2::gsl_specfun::tag for toolbox gsl_specfun
    **/
    struct gsl_sf_gamma_ {};
  }
  NT2_FUNCTION_IMPLEMENTATION(gsl_specfun::tag::gsl_sf_gamma_, gsl_sf_gamma, 1)
  } }

#include <nt2/toolbox/gsl_specfun/functions/scalar/gsl_sf_gamma.hpp>
// #include <nt2/toolbox/gsl_specfun/functions/simd/all/gsl_sf_gamma.hpp> 

#endif

// modified by jt the 29/12/2010

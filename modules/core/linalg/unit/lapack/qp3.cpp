//==============================================================================
//         Copyright 2003 - 2012   LASMEA UMR 6602 CNRS/Univ. Clermont II
//         Copyright 2009 - 2013   LRI    UMR 8623 CNRS/Univ Paris Sud XI
//
//          Distributed under the Boost Software License, Version 1.0.
//                 See accompanying file LICENSE.txt or copy at
//                     http://www.boost.org/LICENSE_1_0.txt
//==============================================================================
#include <nt2/include/functions/qp3.hpp>
#include <nt2/include/functions/fliplr.hpp>
#include <nt2/include/functions/eye.hpp>
#include <nt2/include/functions/zeros.hpp>
#include <nt2/include/functions/colon.hpp>
#include <nt2/include/functions/width.hpp>
#include <nt2/include/functions/height.hpp>
#include <nt2/include/functions/mtimes.hpp>
#include <nt2/include/functions/rand.hpp>

#include <nt2/table.hpp>

#include <nt2/sdk/unit/module.hpp>
#include <nt2/sdk/unit/tests/ulp.hpp>
#include <nt2/sdk/unit/tests/relation.hpp>

NT2_TEST_CASE_TPL(qp3, NT2_REAL_TYPES )
{
  using nt2::_;

  typedef nt2::table<T>           t_t;
  typedef nt2::table<nt2_la_int>  t_i;
  t_t in  = nt2::fliplr( nt2::eye(4, 4, nt2::meta::as_<T>()) );
  t_t tau;
  t_i jpvt;


  t_t ref = -nt2::eye(4, 4, nt2::meta::as_<T>());
  t_t ref_tau = nt2::zeros(4,1 nt2::meta::as_<T>());

  ref_tau(1) = ref_tau(2) = 1;
  ref(4,1)   = ref(3,2)   = 1;

  nt2_la_int p = nt2::qp3(in,jpvt,tau);

  NT2_TEST_EQUAL(in,ref);
  NT2_TEST_EQUAL(tau,ref_tau);


  NT2_TEST_EQUAL(p, 0);
}

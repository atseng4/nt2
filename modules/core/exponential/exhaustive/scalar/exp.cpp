//==============================================================================
//         Copyright 2009 - 2013 LRI    UMR 8623 CNRS/Univ Paris Sud XI
//         Copyright 2013        MetaScale SAS
//
//          Distributed under the Boost Software License, Version 1.0.
//                 See accompanying file LICENSE.txt or copy at
//                     http://www.boost.org/LICENSE_1_0.txt
//==============================================================================
#include <nt2/sdk/unit/exhaustive.hpp>
#include <nt2/include/functions/exp.hpp>

#include <nt2/include/constants/valmin.hpp>
#include <nt2/include/constants/valmax.hpp>
#include <nt2/include/functions/ulpdist.hpp>
#include <nt2/include/functions/mantissa.hpp>
#include <nt2/include/functions/exponent.hpp>

struct raw_exp
{
  float operator()(float x) const
  {
    return (float)nt2::exp(double(x));
  }
};

int main(int ac, char* av[])
{
  float mini = nt2::Valmin<float>();
  float maxi = nt2::Valmax<float>();
  if(ac == 3)
  {
<<<<<<< Updated upstream
    mini = std::atof(av[1]);
    maxi = std::atof(av[2]);
  }
  nt2::exhaustive_test<float> (mini
                              , maxi
                              , nt2::functor<nt2::tag::exp_>()
                              , raw_exp()
                              );
  return 0;
}

################################################################################
#         Copyright 2003 - 2012 LASMEA UMR 6602 CNRS/Univ. Clermont II
#         Copyright 2009 - 2012 LRI    UMR 8623 CNRS/Univ Paris Sud XI
#         Copyright 2012 - 2013 MetaScale SAS
#
#          Distributed under the Boost Software License, Version 1.0.
#                 See accompanying file LICENSE.txt or copy at
#                     http://www.boost.org/LICENSE_1_0.txt
################################################################################
cmake_minimum_required(VERSION 2.8)

set(NT2_DIR ${CMAKE_SOURCE_DIR}/../../../cmake)

find_package(NT2 REQUIRED)
include(${NT2_USE_FILE})

set(CMAKE_BUILD_TYPE Release)
if(CMAKE_COMPILER_IS_GNUCXX)
  add_definitions(-DBOOST_SIMD_NO_STRICT_ALIASING -fno-strict-aliasing)
endif()

add_executable(mandelbrot_nt2 mandelbrot_nt2.cpp)

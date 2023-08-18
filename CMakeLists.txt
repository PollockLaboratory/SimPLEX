cmake_minimum_required( VERSION 3.10 )
set(CMAKE_CXX_STANDARD 17)

project(SimPLEX
  VERSION 0.0
  DESCRIPTION "Bayesian tool for evolutionary analysis of complex substitution patterns"
  LANGUAGES CXX)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall -Wextra")
set(CMAKE_CXX_FLAGS_RELEASE "-w -O3")

add_subdirectory( src )

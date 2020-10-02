cmake_minimum_required( VERSION 3.10 )
set(CMAKE_CXX_STANDARD 17)

project( simplex0_0
  VERSION 0.0
  DESCRIPTION "Tool for evolutionary analysis."
  LANGUAGES CXX)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall -Wextra")
set(CMAKE_CXX_FLAGS_RELEASE "-w -O3")

add_subdirectory( src )

# TESTS
# By default tests are off.
option( BUILD_TESTS "Build the tests." OFF )

if(BUILD_TESTS)
  message(STATUS "Building the tests.")
  # Hard coded path to gtest location. not sure if this is installed with libgtest-dev or googletest package.
  add_subdirectory("/usr/src/gtest" ${CMAKE_BINARY_DIR}/gtest)
  enable_testing()
  include_directories(${GTEST_INCLUDE_DIRS})
  add_subdirectory(tests)
endif()


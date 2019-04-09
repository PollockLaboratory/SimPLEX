cmake_minimum_required( VERSION 3.0 )
set(CMAKE_CXX_STANDARD 17)

project( simplex0_0 )

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

find_package(Boost 1.56 REQUIRED)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
endif()

find_package(Lua 5.2 REQUIRED)
include_directories(${LUA_INCLUDE_DIR})

add_subdirectory( src )
include_directories( src libs )


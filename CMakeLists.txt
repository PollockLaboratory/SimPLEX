cmake_minimum_required( VERSION 3.0 )
set(CMAKE_CXX_STANDARD 17)
project( simplex0_0 )

find_package(Lua 5.2 REQUIRED)
include_directories(${LUA_INCLUDE_DIR})

add_subdirectory( src )
include_directories( src libs )


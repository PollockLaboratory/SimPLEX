# Add sources.
set(SRCS
  Data.cpp
  MCMC.cpp
  Environment.cpp
  IO/Files.cpp
  IO/TreeParser.cpp
  IO/SubstitutionModelParser.cpp
  IO/SequencesParser.cpp
  IO/SubstitutionModelParameterWrapper.cpp
  IO/LuaUtils.cpp
  SubstitutionCounts.cpp
  Model.cpp
  ModelParts/Sequence.cpp
  ModelParts/ComponentSet.cpp
  ModelParts/AbstractComponent.cpp
  ModelParts/SubstitutionModels/SubstitutionModel.cpp
  ModelParts/SubstitutionModels/RateVector.cpp
  ModelParts/SubstitutionModels/Parameters.cpp
  ModelParts/Trees/Tree.cpp
  ModelParts/Trees/BranchSplitting.cpp
  ModelParts/Trees/TreeParts.cpp )

add_library(libSimPLEX STATIC ${SRCS})

find_package(Boost 1.5 REQUIRED)
if(Boost_FOUND)
  target_include_directories( libSimPLEX PUBLIC ${Boost_INCLUDE_DIRS})
  target_link_libraries( libSimPLEX PUBLIC ${Boost_LIBRARIES})
endif()

find_package(Lua 5 REQUIRED)
if(Lua_FOUND)
  target_include_directories( libSimPLEX PUBLIC ${LUA_INCLUDE_DIR})
  target_link_libraries( libSimPLEX PUBLIC ${LUA_LIBRARIES} )
endif()

target_include_directories( libSimPLEX PUBLIC ../libs )

# Main executable.
add_executable( simPLEX SimPLEX.cpp)
target_link_libraries( simPLEX PUBLIC libSimPLEX )
set_target_properties( simPLEX PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin" )

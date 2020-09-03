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

add_executable( simPLEX SimPLEX.cpp)

find_package(Boost 1.5 REQUIRED)
if(Boost_FOUND)
  target_include_directories( simPLEX PUBLIC ${Boost_INCLUDE_DIRS})
  target_link_libraries( simPLEX PUBLIC ${Boost_LIBRARIES})
endif()

find_package(Lua 5 REQUIRED)
if(Lua_FOUND)
  target_include_directories( simPLEX PUBLIC ${LUA_INCLUDE_DIR})
  target_link_libraries( simPLEX PUBLIC ${LUA_LIBRARIES} )
endif()

target_sources( simPLEX PUBLIC ${SRCS} )

target_include_directories( simPLEX PUBLIC ../libs )

set_target_properties( simPLEX PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin" )

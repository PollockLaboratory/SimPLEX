add_subdirectory( SubstitutionModels )
include_directories( SubstitutionModels )

set(SRCS
  utils.cpp
  Data.cpp
  MCMC.cpp
  Environment.cpp
  IO/Files.cpp
  IO/TreeParser.cpp
  IO/SubstitutionModelParser.cpp
  IO/SequencesParser.cpp
  Sequence.cpp
  SubstitutionCounts.cpp
  Model.cpp
  SubstitutionModels/SubstitutionModel.cpp
  SubstitutionModels/Components/AbstractComponent.cpp
  SubstitutionModels/Components/ComponentSet.cpp
  SubstitutionModels/Components/ComponentTypes.cpp
  SubstitutionModels/Components/RateVector.cpp
  SubstitutionModels/Components/AbstractValueTypes.cpp
  SubstitutionModels/Components/SampleableValueTypes.cpp
  Trees/Tree.cpp
  Trees/BranchSplitting.cpp
  Trees/TreeParts.cpp )

add_executable( simPLEX SimPLEX.cpp)
target_sources( simPLEX PUBLIC ${SRCS} )
target_link_libraries( simPLEX ${LUA_LIBRARIES} )
target_link_libraries( simPLEX ${Boost_LIBRARIES})

include_directories( ../libs )

message( STATUS "${LUA_LIBRARIES}" )
message( STATUS "${LUA_INCLUDE_DIR}" )

set_target_properties( simPLEX PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin" )

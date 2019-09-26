set(SRCS
  Data.cpp
  MCMC.cpp
  Environment.cpp
  IO/Files.cpp
  IO/TreeParser.cpp
  IO/SubstitutionModelParser.cpp
  IO/SequencesParser.cpp
  SubstitutionCounts.cpp
  Model.cpp
  ModelParts/Sequence.cpp
  ModelParts/ComponentSet.cpp
  ModelParts/AbstractComponent.cpp
  ModelParts/SubstitutionModels/SubstitutionModel.cpp
  ModelParts/SubstitutionModels/Components/ComponentTypes.cpp
  ModelParts/SubstitutionModels/Components/RateVector.cpp
  ModelParts/SubstitutionModels/Components/AbstractValueTypes.cpp
  ModelParts/SubstitutionModels/Components/SampleableValueTypes.cpp
  ModelParts/Trees/Tree.cpp
  ModelParts/Trees/BranchSplitting.cpp
  ModelParts/Trees/TreeParts.cpp )

add_executable( simPLEX SimPLEX.cpp)
target_sources( simPLEX PUBLIC ${SRCS} )
target_link_libraries( simPLEX ${LUA_LIBRARIES} )
target_link_libraries( simPLEX ${Boost_LIBRARIES})

include_directories( ../libs )

message( STATUS "${LUA_LIBRARIES}" )
message( STATUS "${LUA_INCLUDE_DIR}" )

set_target_properties( simPLEX PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin" )

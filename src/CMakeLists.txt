add_subdirectory( SubstitutionModels )
include_directories( SubstitutionModels )

set(SRCS
  utils.cpp
  Data.cpp
  MCMC.cpp
  Environment.cpp
  IO.cpp
  Sequence.cpp
  SubstitutionCounts.cpp
  Model.cpp
  SubstitutionModels/SubstitutionModel.cpp
  SubstitutionModels/SubstitutionModelTypes.cpp
  SubstitutionModels/Components/AbstractComponent.cpp
  SubstitutionModels/Components/ComponentSet.cpp
  SubstitutionModels/Components/ComponentTypes.cpp
  SubstitutionModels/Components/RateVector.cpp
  SubstitutionModels/Components/Types/AbstractValueTypes.cpp
  SubstitutionModels/Components/Types/SampleableValueTypes.cpp
  SubstitutionModels/Types/CategoryGTR.cpp
  SubstitutionModels/Types/FixedRate.cpp
  SubstitutionModels/Types/GeneralTimeReversible.cpp
  SubstitutionModels/Types/SingleRate.cpp
  Trees/Tree.cpp
  Trees/TreeParts.cpp
  Trees/BranchSplitting.cpp
  Trees/TreeParser.cpp )

add_executable( simPLEX SimPLEX.cpp)
target_sources( simPLEX PUBLIC ${SRCS} )

set_target_properties( simPLEX PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin" )

# This holds all the sources that I want to include into SimPLEX
# It is here so I can have other .cpp files in the directory without them 
# being compiled and linked in. For example, I have a main.cpp that is a 
# separate testing file and is not part of SimPLEX. 
 
SOURCES = \
 SimPLEX.cpp \
 Options.cpp \
 Data.cpp \
 MCMC.cpp \
 Model.cpp \
 Tree.cpp \
 Tree_B1.cpp \
 SubstitutionModel.cpp \
 SingleProbabilityModel.cpp \
 SingleSubstitutionRateModel.cpp \
 ProbabilityMatrixModel.cpp \
 SubstitutionRateMatrixModel.cpp \
 RadiusDependentModel.cpp \
 SubstitutionMixtureModel.cpp \
 
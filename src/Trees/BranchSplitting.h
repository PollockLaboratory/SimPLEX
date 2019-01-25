#ifndef BranchSplitting_h_
#define BranchSplitting_h_

#include <utility>
#include <functional>
#include "TreeParts.h"
#include <cmath> // for floor and pow

// Dispatch function - choose the correct algorithm to return to Tree class.
std::function< std::pair<BranchSegment*, BranchSegment*>(float)> pickBranchSplitAlgorithm();

// Possible algorithms.
std::pair<BranchSegment*, BranchSegment*> noSplitMethod(float distance); 
std::pair<BranchSegment*, BranchSegment*> splitHalfMethod(float distance); 

#endif



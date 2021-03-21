#include "BranchSplitting.h"
#include "../../Environment.h"

extern Environment env;

// Dispatch function.
std::function< std::pair<BranchSegment*, BranchSegment*>(float)> pickBranchSplitAlgorithm() {
  int option = env.get<int>("TREE.branch_split_algorithm");
  std::function< std::pair<BranchSegment*, BranchSegment*>(float)> f;
  if(option == 0) {
    f = noSplitMethod;
  } else if(option == 1) {
    f = splitHalfMethod;
  } else {
    std::cout << "Error: Invalid option for branch_split_algorithn: " << option << " Option not recognized." << std::endl;
    exit(EXIT_FAILURE);
  }
  return(f);
}

// Possible algorithms.
std::pair<BranchSegment*, BranchSegment*> noSplitMethod(float distance) {
  BranchSegment* newBranchTop = new BranchSegment(distance);
  BranchSegment* newBranchBottom = newBranchTop;

  return(std::make_pair(newBranchTop, newBranchBottom));
}

std::pair<BranchSegment*, BranchSegment*> splitHalfMethod(float distance) {
  /* 
   * This algorithm works by splitting branches in half until each segment is less
   * than the max sequence length.
   */

  static double u = env.get<double>("UNIFORMIZATION.initial_value");
  static double max_segment_probability = env.get<double>("TREE.max_segment_probability");

  //static double max_seg_len = max_segment_probability / ((1.0 - max_segment_probability) * u);
  static double max_seg_len = max_segment_probability / u;
  
  float dist = distance;

  // Calculate how many internal nodes are needed for given branch.
  int extraNodes = 0;
  for(int i = 0; dist > max_seg_len; i++) {
    dist = dist/2.0;
    extraNodes += pow(2, i);
  }

  BranchSegment* newBranchTop = new BranchSegment(dist);
  BranchSegment* newBranchBottom = newBranchTop;

  //branchList.push_back(newBranchTop); //Adds the new branch to the branchList.

  //Create extra nodes and link together.
  if(extraNodes > 0) {
    for(int i = 0; i < extraNodes; i++) {
      TreeNode* n = new TreeNode();
      n->distance = dist;

      BranchSegment* b = new BranchSegment(dist);

      newBranchBottom->decendant = n;
      n->up = newBranchBottom;

      n->left = b;
      b->ancestral = n;

      newBranchBottom = b;
    }
  }
  return(std::make_pair(newBranchTop, newBranchBottom));
}



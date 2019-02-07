#ifndef TreeParts_h_
#define TreeParts_h_

#include <string>
#include <map>
#include <vector>
#include <iostream>
#include <math.h>

#include "TreeParser.h"
#include "../Sequence.h"
#include "Components/RateVector.h"
#include "../SubstitutionModels/SubstitutionModel.h"

class TreeNode;

class BranchSegment {
 public:
  BranchSegment(float distance);
  ~BranchSegment();

  float distance;
  TreeNode* ancestral;
  TreeNode* decendant;

  double get_rate(int pos, int dec_state);

  // Key statistics.
  inline void update_rate_vectors();
  bool virtualSubstituionQ(int state);

  void update_counts(std::map<RateVector*, std::vector<int>>& subs_by_rateVector, std::pair<int, int>& subs_by_branch);

  void update();
 
  friend std::ostream& operator<< (std::ostream &out, const BranchSegment &b);
 private:
  std::vector<RateVector*> rates; //By site.
};

class TreeNode {
 public:
  static int unique_id;
  std::string name;
  double distance;

  BranchSegment* up;
  BranchSegment* left;
  BranchSegment* right;

  std::vector<int>* sequence;
  SequenceAlignment* MSA;

  SubstitutionModel* SM;

  TreeNode();
  TreeNode(IO::RawTreeNode* raw_tree);
  TreeNode(std::string n);

  void sample();
  void sampleSequence();
  void sampleSinglePosition(int pos);
  std::string toString();

  bool isTip();
  bool sampled;
};

#endif

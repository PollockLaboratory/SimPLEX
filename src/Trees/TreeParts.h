#ifndef TreeParts_h_
#define TreeParts_h_

#include <string>
#include <map>
#include <vector>
#include <iostream>
#include <math.h>

#include "../IO/TreeParser.h"
#include "../Sequence.h"
#include "Components/RateVector.h"
#include "../SubstitutionModels/SubstitutionModel.h"
#include "../SubstitutionCounts.h"

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
  void set_new_substitutions();

  void update_counts(std::map<RateVector*, std::vector<int>>& subs_by_rateVector,
		       raw_counts& subs_by_branch);

  void update();
 
  friend std::ostream& operator<< (std::ostream &out, const BranchSegment &b);

  std::vector<bool> substitutions;
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

  TreeNode* sample();
  void sample_sequence();
  int sample_single_position(int pos);

  // Printing/Display
  std::string toString();
  std::string get_sequence();
  std::string state_at_pos(int i);

  bool isTip();
  bool sampled;
};

#endif

#ifndef Tree_h_
#define Tree_h_

#include <algorithm>
#include <functional>
#include <istream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <list>

#include "Sequence.h"
#include "SubstitutionModel.h"
#include "Trees/TreeParser.h"
#include "TreeParts.h"
#include "BranchSplitting.h"
#include "SubstitutionCounts.h"

using std::string;
using std::map;
using std::vector;

class Tree {
 public:
  TreeNode* root;
  SequenceAlignment* MSA;
  SubstitutionModel* SM;

  map<string, vector<int>> names_to_sequences;

  std::list<BranchSegment*> branchList; // Potentially should be vectors.
  std::vector<TreeNode*> nodeList;

  // Constructing/Initializing.
  Tree();
  void Initialize(IO::RawTreeNode* raw_tree, SequenceAlignment* &MSA, SubstitutionModel* &SM);

  // Internal Tree nodes.
  void connect_nodes(TreeNode* &ancestral, BranchSegment* &ancestralBP, TreeNode* &decendant, float distance);
  TreeNode* createTreeNode(IO::RawTreeNode* raw_tree, TreeNode* &ancestralNode, BranchSegment* &ancestralBP);

  // Sampling.
  bool sample();

  // Possible sampling methods.
  bool(Tree::*treeSamplingMethod)(); 
  bool sample_ancestral_states(); // When the tree is actually being sampled.
  bool step_through_MSAs(); // When the ancestral sequences have already been determined. 

  // Counts
  void update_counts(SubstitutionCounts&); // New.
  std::list<float> get_branch_lengths();

  // Recording state data.
  static std::ofstream tree_out;
  static std::ofstream substitutions_out;
  void initialize_output_streams();
  void record_tree();
  void record_state(int gen, double l);

  // Debug tools.
  void print_branchList();
  void print_nodeList();
  void print_parameters();
 private:
  // Settings/options.
  float max_seg_len; // Max segment length.
  int seqLen;
  std::function< std::pair<BranchSegment*, BranchSegment*>(float)> splitBranchMethod; // Algorithm for splitting branches.
  float u;

  void configureSequences(TreeNode* n);
};

#endif

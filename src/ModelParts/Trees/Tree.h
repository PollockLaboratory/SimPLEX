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

#include "../Sequence.h"
#include "../SubstitutionModels/SubstitutionModel.h"
#include "../../IO/TreeParser.h"
#include "TreeParts.h"
#include "BranchSplitting.h"
#include "../../SubstitutionCounts.h"

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
  std::list<TreeNode*> tipList;

  // Constructing/Initializing.
  Tree();
  void Initialize(IO::RawTreeNode* raw_tree, SequenceAlignment* &MSA, SubstitutionModel* &SM);

  // Internal Tree nodes.
  void connect_nodes(TreeNode* &ancestral, BranchSegment* &ancestralBP, TreeNode* &decendant, float distance);
  TreeNode* createTreeNode(IO::RawTreeNode* raw_tree, TreeNode* &ancestralNode, BranchSegment* &ancestralBP);

  // Mark gaps.
  void identify_gaps();
  
  // Sampling.
  sample_status sample();
  void find_ancestral_sequences();

  // Possible sampling methods.
  sample_status(Tree::*treeSamplingMethod)(); 
  sample_status sample_ancestral_states(); // When the tree is actually being sampled.
  sample_status step_through_MSAs(); // When the ancestral sequences have already been determined. 

  // Counts
  void update_counts(SubstitutionCounts&) const; // New.
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
private:
  // Settings/options.
  float max_seg_len; // Max segment length.
  int seqLen;
  std::function< std::pair<BranchSegment*, BranchSegment*>(float)> splitBranchMethod; // Algorithm for splitting branches.

  void configureSequences(TreeNode* n);

  // Recording State.
  void record_substitutions(int gen, double l);
};

class TreeParameter : public SampleableComponent {
private:
  Tree* tree;
public:
  TreeParameter();

  virtual void print();
  
  virtual sample_status sample();
  virtual void undo();
  virtual void fix();
  virtual void refresh();
  virtual double record_state(int gen, double l);

  void Initialize(IO::RawTreeNode* raw_tree, IO::RawMSA* &raw_msa, SubstitutionModel* &SM);
  Tree* get_tree_ptr();
};

#endif

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
private:
  TreeNode* root;
  SequenceAlignment* MSA;
  SubstitutionModel* SM;
  map<string, vector<int>> names_to_sequences;

  std::list<BranchSegment*> branchList; // Potentially should be vectors.
  std::vector<TreeNode*> nodeList;
  std::list<TreeNode*> tipList;

  // Settings/options.
  float max_seg_len; // Max segment length.
  int seqLen;
  std::function< std::pair<BranchSegment*, BranchSegment*>(float)> splitBranchMethod; // Algorithm for splitting branches.

  // Initialize.
  void connect_nodes(TreeNode* &ancestral, BranchSegment* &ancestralBP,
		     TreeNode* &decendant, float distance);
  TreeNode* createTreeNode(IO::RawTreeNode* raw_tree, TreeNode* &ancestralNode,
			   BranchSegment* &ancestralBP);
  void identify_gaps();
  void configureSequences(TreeNode* n);

  // Output.
  void record_substitutions(int gen, double l);

  // Possible sampling methods.
  sample_status(Tree::*treeSamplingMethod)(const std::list<int>&); 
  sample_status sample_ancestral_states(const std::list<int>&); // When the tree is actually being sampled.
  sample_status step_through_MSAs(const std::list<int>&); // When the ancestral sequences have already been determined. 
public:
  // Constructing/Initializing.
  Tree();
  void Initialize(IO::RawTreeNode* raw_tree, SequenceAlignment* &MSA, SubstitutionModel* &SM);
  sample_status sample(const std::list<int>& positions);
  
  SubstitutionModel* get_SM();
  const std::list<BranchSegment*> get_branches();
  std::list<float> get_branch_lengths();

  void record_tree();
  void record_state(int gen, double l);

  // Debug tools.
  void print_branchList();
  void print_nodeList();
};

class AncestralStatesParameter : public SampleableComponent {
private:
  Tree* tree;
  int n_samples; // Number of positions/columns to resample each sample request.
public:
  AncestralStatesParameter();

  void print() override;
  std::string get_type() override;
  
  sample_status sample() override;
  void undo() override;
  void fix() override;
  void refresh() override;
  double record_state(int gen, double l) override;

  void Initialize(IO::RawTreeNode* raw_tree, IO::RawMSA* &raw_msa, SubstitutionModel* &SM);
  Tree* get_tree_ptr();
};

#endif

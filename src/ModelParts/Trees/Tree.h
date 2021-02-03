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
  SubstitutionModel* SM;
  map<string, vector<int>> names_to_sequences;
  
  // Alternative node access.
  void buildNodeLists(TreeNode*);
  std::list<BranchSegment*> branchList; // Potentially should be vectors.
  std::list<TreeNode*> nodeList;
  std::list<TreeNode*> tipList;

  // Settings/options.
  std::function< std::pair<BranchSegment*, BranchSegment*>(float)> splitBranchMethod; // Algorithm for splitting branches.

  // Initialize.
  void connect_nodes(TreeNode* &ancestral, BranchSegment* &ancestralBP,
		     TreeNode* &decendant, float distance);
  TreeNode* createTreeNode(IO::RawTreeNode* raw_tree, TreeNode* &ancestralNode,
			   BranchSegment* &ancestralBP, float scale_factor);

  // Possible sampling methods.
  sample_status(Tree::*treeSamplingMethod)(const std::list<int>&); 
  sample_status sample_ancestral_states(const std::list<int>&); // When the tree is actually being sampled.
  sample_status step_through_MSAs(const std::list<int>&); // When the ancestral sequences have already been determined. 
public:
  // Constructing/Initializing.
  Tree();
  void Initialize(IO::RawTreeNode* raw_tree);
  sample_status sample(const std::list<int>& positions);
  
  SubstitutionModel* get_SM();
  const std::list<BranchSegment*> get_branches();
  const std::list<TreeNode*> nodes();
  std::list<float> get_branch_lengths();

  void record_tree();
  void record_state(int gen, double l);

  // Output.
  void record_substitutions(int gen, double l);

  // Debug tools.
  void print_branchList();
  void print_nodeList();

  // Sequence stuff.
  TreeNode* root;
  SequenceAlignment* MSA;
  void configureSequences(TreeNode* n);
  void configureBranches(TreeNode* n, unsigned int n_columnes);
  void connect_substitution_model(SubstitutionModel*);
};

class RateVectorAssignmentParameter : public AbstractComponent {
private:
  Tree* tree;
public:
  RateVectorAssignmentParameter(Tree*);

  void print() override;
  std::string get_type() override;
  
  void fix() override;
  void refresh() override;

  std::string get_state_header() override;
  std::string get_state() override; 
};

#endif

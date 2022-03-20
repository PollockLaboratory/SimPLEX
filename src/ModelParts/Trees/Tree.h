#ifndef Tree_h_
#define Tree_h_

#include <functional>
#include <string>
#include <map>
#include <vector>
#include <list>

#include "../Sequence.h"
#include "../SubstitutionModels/SubstitutionModel.h"
#include "../../IO/TreeParser.h"
#include "TreeParts.h"
#include "BranchSplitting.h"

using std::string;
using std::map;
using std::vector;

class Tree {
private:
  map<string, vector<int>> names_to_sequences;
  
  // Cache fast node access.
  void cache_node_access();
  void buildNodeLists(TreeNode*);
  void build_recursion_path(TreeNode*, std::list<TreeNode*>&);

  std::list<BranchSegment*> branchList; // Potentially should be vectors.
  std::list<TreeNode*> nodeList; // In order of tree recursion.
  std::list<TreeNode*> tipList;

  std::vector<TreeNode*> nodeVector;

  std::map<TreeNode*, std::list<TreeNode*>> recursion_paths;

  // Settings/options.
  std::function< std::pair<BranchSegment*, BranchSegment*>(float)> splitBranchMethod; // Algorithm for splitting branches.

  // Initialize.
  void connect_nodes(TreeNode* &ancestral, BranchSegment* &ancestralBP,
		     TreeNode* &decendant, float distance);
  TreeNode* createTreeNode(IO::RawTreeNode* raw_tree, TreeNode* &ancestralNode,
			   BranchSegment* &ancestralBP, float scale_factor);
  
public:
  SubstitutionModel* SM; // Temp - should be private.

  // Constructing/Initializing.
  Tree();
  void Initialize(IO::RawTreeNode* raw_tree);
  
  SubstitutionModel* get_SM();

  const std::list<BranchSegment*>& get_branches();
  const std::list<TreeNode*>& nodes();
  const std::list<TreeNode*>& get_recursion_path(TreeNode*);

  std::list<float> get_branch_lengths();

  TreeNode* rand_node();

  //Output
  void record_tree();

  // Debug tools.
  void print_branchList();
  void print_nodeList();

  // Sequence stuff.
  TreeNode* root;
  void configure_branches(TreeNode* n, unsigned int n_columnes, std::list<std::string> state_domain_names);
  void connect_substitution_model(SubstitutionModel*);
};

// Parameter that triggers a refresh of all the RateVector assignments on the Branch segments.
// Also updates the substitutions on the branch segments as well.
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

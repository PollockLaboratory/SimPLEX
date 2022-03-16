#include "Tree.h"

#include <iostream>
#include <sstream> // For ostringstream
#include <cmath> // for floor and pow
#include <unordered_set>

#include "../../Environment.h"
#include "../../IO/Files.h"

extern double Random();
extern Environment env;
extern IO::Files files;

// Tree constructor.
Tree::Tree() {   
}

// Tree Initialize using seqs and states
void Tree::Initialize(IO::RawTreeNode* raw_tree) {
  // Set dynamicly chosen branch-splitting function.
  splitBranchMethod = pickBranchSplitAlgorithm();

  // Scaling the tree.
  float scale_factor = 1.0;
  if(env.get<bool>("TREE.scale_tree")){ 
    float total_tree_length = IO::findRawTreeTotalLength(raw_tree);
    float target_length = env.get<double>("TREE.target_tree_length");
    scale_factor = target_length/total_tree_length;
    std::cout << "\tRescaling tree from total length " << total_tree_length << " to " << target_length << std::endl;
  } else {
    scale_factor = 1.0;
  }
 
  //this->MSA = MSA;
  //this->SM = SM;
  //SM->organizeRateVectors(MSA->numCols(), SM->get_states()->n);

  // Proxys are created to correctly create root node.
  BranchSegment* proxyBranch = new BranchSegment(0.0);
  TreeNode* proxyNode = new TreeNode("Proxy");
  root = createTreeNode(raw_tree, proxyNode, proxyBranch, scale_factor);
  delete proxyBranch;
  delete proxyNode;

  buildNodeLists(root);

  std::cout << "\tTree split into " << nodeList.size() << " nodes." << std::endl;

  //Setup output.
  files.add_file("tree_out", env.get<std::string>("OUTPUT.tree_out_file"), IOtype::OUTPUT);
  record_tree();

  //files.add_file("substitutions_out", env.get<std::string>("OUTPUT.substitutions_out_file"), IOtype::OUTPUT);
  //files.write_to_file("substitutions_out", "I,GEN,LogL,Ancestral,Decendant,Substitutions\n");

}

// Creation of tree nodes.
void Tree::connect_nodes(TreeNode* &ancestralNode, BranchSegment* &ancestralBP, TreeNode* &decendantNode, float distance) {
  /* This function is responsible for connecting two Nodes with branch segments.
   * This includes braking up long branches into smaller branch segment when needed and creating new internal branch nodes.
   * Arguments:
   *    ancestralNode - pointer to the ancestral node.
   * 	ancestralBP - pointer of the ancestral node that points to the decendant node. (This is needed to distinguish between the left and right branch pointers of ancestral node.)
   * 	decendantNode - pointer to the decendantNode.
   * 	distance - the branch length.
   */
  std::pair<BranchSegment*, BranchSegment*> intermediateBranches = splitBranchMethod(distance);
  BranchSegment* newBranchTop = intermediateBranches.first;
  BranchSegment* newBranchBottom = intermediateBranches.second;

  decendantNode->up = newBranchBottom;
  ancestralBP = newBranchTop;

  newBranchTop->ancestral = ancestralNode;
  newBranchBottom->decendant = decendantNode;
}

TreeNode* Tree::createTreeNode(IO::RawTreeNode* raw_tree, TreeNode* &ancestralNode, BranchSegment* &ancestralBP, float scale_factor) {
  TreeNode* newTreeNode = new TreeNode(raw_tree);

  connect_nodes(ancestralNode, ancestralBP, newTreeNode, raw_tree->distance * scale_factor);
  newTreeNode->distance = ancestralBP->distance; // Correct tree node distance for splitting.

  if(raw_tree->left != 0) {
    createTreeNode(raw_tree->left, newTreeNode, newTreeNode->left, scale_factor);
  }

  if(raw_tree->right != 0) {
    createTreeNode(raw_tree->right, newTreeNode, newTreeNode->right, scale_factor);
  }

  return(newTreeNode);
}

void Tree::buildNodeLists(TreeNode* n) {
  nodeList.push_front(n);

  if(n->left != 0) {
    BranchSegment* b = n->left;
    branchList.push_back(b);
    buildNodeLists(b->decendant);
  }

  if(n->right != 0) {
    BranchSegment* b = n->right;
    branchList.push_back(b);
    buildNodeLists(b->decendant);
  }

  if(n->isTip()) {
      tipList.push_back(n);
  }
}

void Tree::configureBranches(TreeNode* n, unsigned int n_columns, std::map<std::string, std::list<std::string>> all_states) {
  if(n->left != 0) {
    BranchSegment* b = n->left;
    b->Initialize(n_columns, all_states);
    configureBranches(b->decendant, n_columns, all_states);
  }

  if(n->right != 0) {
    BranchSegment* b = n->right;
    b->Initialize(n_columns, all_states);
    configureBranches(b->decendant, n_columns, all_states);
  }
}

void Tree::connect_substitution_model(SubstitutionModel* sm) {
  this->SM = sm;

  for(auto n = nodeList.begin(); n != nodeList.end(); ++n) {
    (*n)->connect_substitution_model(sm);
  }
}

SubstitutionModel* Tree::get_SM() {
  return(SM);
}

// Sampling and likelihood.
const std::list<BranchSegment*> Tree::get_branches() {
  return(branchList);
}

const std::list<TreeNode*> Tree::nodes() {
  return(nodeList);
}

std::list<float> Tree::get_branch_lengths() {
  // Maybe Tree should just hold onto all the branch lengths in play?
  std::list<float> lens = {};
  std::unordered_set<float> lens_set = {};
  BranchSegment* b;
  for(auto it = branchList.begin(); it != branchList.end(); ++it) {
    b = *it;
    if(lens_set.find(b->distance) == lens_set.end()) {
      // Branch does NOT already exists.
      lens_set.insert(b->distance);
      lens.push_back(b->distance);
    }
  }
  return(lens);
}

TreeNode* Tree::rand_node() {
  // NOTE probably a faster way to do this, but its not a critical function call.
  auto it = nodeList.begin();
  for(int offset = rand() % nodeList.size(); offset > 0; offset--) {
    it++;
  }
  return(*it);
}

// Record State data.
void Tree::record_tree() {
  /*
   * Records the tree topology with all node names and branch segments.
   */
  std::string tree_str = root->toString() + ";";
  files.write_to_file("tree_out", tree_str);
}

// Debug tools.
void Tree::print_branchList() {
  std::cout << "Printing Branch list. Size: " << branchList.size() << std::endl;
  for(std::list<BranchSegment*>::iterator it = branchList.begin(); it != branchList.end(); ++it) {
    BranchSegment* b = *it;
    std::cout << "Branch: " << b->distance << std::endl;
  }
}

void Tree::print_nodeList() {
  std::cout << "Printing Node list. Size:  " << nodeList.size() << std::endl;
  for(auto it = nodeList.begin(); it != nodeList.end(); ++it) {
    std::cout << "Node: " << (*it)->name << std::endl;
  }
}

// RateVectorAssignmentParameter

RateVectorAssignmentParameter::RateVectorAssignmentParameter(Tree* tree) : AbstractComponent("RateVectorassignment") {
  this->tree = tree;
}

void RateVectorAssignmentParameter::print() {
  std::cout << "RateVectorAssignment" << std::endl;
}

std::string RateVectorAssignmentParameter::get_type() {
  return("RateVectorAssignment");
}
 
void RateVectorAssignmentParameter::fix() {
}

void RateVectorAssignmentParameter::refresh() {
  // Update all branches - new substitutions.
  std::list<BranchSegment*> branches = tree->get_branches();
  for(auto b = branches.begin(); b != branches.end(); ++b) {
    (*b)->update();
  }
}

std::string RateVectorAssignmentParameter::get_state_header() {
  return(name);
}

std::string RateVectorAssignmentParameter::get_state() {
  return("n/a");
}


#include "Tree.h"

#include <iostream>
#include <sstream> // For ostringstream
#include <cmath> // for floor and pow
#include <unordered_set>

#include "../Sequence.h"
#include "../Environment.h"
#include "../IO.h"
#include "../utils.h"

extern double Random();
extern Environment env;
extern IO::Files files;

using namespace std;

ofstream Tree::substitutions_out;
ofstream Tree::tree_out;

// Tree constructor.
Tree::Tree() {   
}

// Tree Initialize using seqs and states
void Tree::Initialize(IO::RawTreeNode* raw_tree, SequenceAlignment* &MSA, SubstitutionModel* &SM) {
  std::cout << "Creating MCMC tree structure." << std::endl;

  // Configuration
  max_seg_len = env.get_float("max_segment_length");
  u = env.get_float("uniformization_constant");

  // Dynamicly chosen functions.
  splitBranchMethod = pickBranchSplitAlgorithm();

  if(env.ancestral_sequences) {
    treeSamplingMethod = &Tree::step_through_MSAs;
  } else {
    treeSamplingMethod = &Tree::sample_ancestral_states;
  }
  
  this->MSA = MSA;
  seqLen = MSA->numCols();
  this->SM = SM;

  // Proxys are created to correctly create root node.
  BranchSegment* proxyBranch = new BranchSegment(0.0);
  TreeNode* proxyNode = new TreeNode("Proxy");
  root = createTreeNode(raw_tree, proxyNode, proxyBranch);
  delete proxyBranch;
  delete proxyNode;

  std::cout << "Attaching sequences to tree." << std::endl;
  configureSequences(root);
  
  for(auto t = nodeList.begin(); t != nodeList.end(); ++t) {
    if((*t)->isTip()) {
      tipList.push_back(*t);
    }
  }

  // Initial sample to get counts.
  sample_ancestral_states();

  for(auto b = branchList.begin(); b != branchList.end(); ++b) {
  	(*b)->update();
  }

  initialize_output_streams();	

  std::cout << std::endl;
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

TreeNode* Tree::createTreeNode(IO::RawTreeNode* raw_tree, TreeNode* &ancestralNode, BranchSegment* &ancestralBP) {
  TreeNode* newTreeNode = new TreeNode(raw_tree);
  connect_nodes(ancestralNode, ancestralBP, newTreeNode, raw_tree->distance);
  newTreeNode->distance = ancestralBP->distance; // Correct tree node distance for splitting.

  if(raw_tree->left != 0) {
    createTreeNode(raw_tree->left, newTreeNode, newTreeNode->left);
  }

  if(raw_tree->right != 0) {
    createTreeNode(raw_tree->right, newTreeNode, newTreeNode->right);
  }

  return(newTreeNode);
}

// Configure sequences.
void Tree::configureSequences(TreeNode* n) {
  /*
   * Traverses the tree attaching sequences to nodes.
   * Also adds all the Nodes and Branch segments to their corresponding lists.
   */
  n->MSA = MSA;
  n->SM = SM;

  // Add nodes to nodeList and BranchSegments to branchList.
  // Recursively call cofigureSequences on rest of the tree.
  nodeList.push_back(n);

  if(n->left != 0) {
    BranchSegment* b = n->left;
    branchList.push_back(b);
    configureSequences(b->decendant);
  }

  if(n->right != 0) {
    BranchSegment* b = n->right;
    branchList.push_back(b);
    configureSequences(b->decendant);
  }

  if(MSA->taxa_names_to_sequences.count(n->name)) {
    n->sequence = &(MSA->taxa_names_to_sequences.at(n->name));
  } else {
    if(env.ancestral_sequences) {
      // ANCESTRAL SEQUENCES KNOWN
      // In the case when ancestral sequences are known there can be no missing sequences.
      // New sequences should not be created.
      std::cerr << "Error: Missing sequence for \"" << n->name << "\"." << std::endl;
      exit(EXIT_FAILURE);
    } else {
      // NORMAL RUN
      // only tip sequences are needed.
      if(n->isTip()){
	std::cerr << "Error: Missing sequence for \"" << n->name << "\"." << std::endl;
	exit(EXIT_FAILURE);
      } else {
	// Add new sequence to sequence alignments.
	MSA->add(n->name);
	n->sequence = &(MSA->taxa_names_to_sequences.at(n->name));

	// Fill missing sequences/
	if(n->left != 0 and n->right == 0) {
	  // Internal Continous.
	  TreeNode* dsNode = n->left->decendant; // ds = downstream.
	  *(n->sequence) = *(dsNode->sequence);
	} else {
	  // Root or internal branch.
	  vector<int> dsNodeLseq = *(n->left->decendant->sequence);
	  vector<int> dsNodeRseq = *(n->right->decendant->sequence);
	  vector<int> p = MSA->findParsimony(dsNodeLseq, dsNodeRseq);
	  *(n->sequence) = p;
	}
      }
    }
  }
}

// Sampling and likelihood.
void Tree::update_counts(SubstitutionCounts& counts) {
  for(auto it = branchList.begin(); it != branchList.end(); ++it) {
    BranchSegment* b = *it;
    b->update_counts(counts.subs_by_rateVector, counts.subs_by_branch[b->distance]);
  }
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

//
// SAMPLING TREE PARAMETERS.
//

bool Tree::sample() {
  bool s = (this->*treeSamplingMethod)();

  // Update branch list - new substitutions.
  for(auto b = branchList.begin(); b != branchList.end(); ++b) {
    (*b)->update();
  }

  return(s);
}

//
// Two options for changing ancestral states.
//

void Tree::traverse_find_ancestral_sequences() {
  std::queue<TreeNode*> nodes = {};
  for(auto t = tipList.begin(); t != tipList.end(); ++t) {
    (*t)->sampled = true;
    nodes.push((*t)->up->ancestral);
  }

  while(not nodes.empty()) {
    if(not nodes.front()->sampled) {
      nodes.push(nodes.front()->sample());
    }

    nodes.pop();
  }

  for(auto n = nodeList.begin(); n != nodeList.end(); ++n) {
    (*n)->sampled = false;
  }
}

bool Tree::sample_ancestral_states() {
  /*
   * Both recalculates substitutions and recalculates the ancestral sequences.
   */

  for(auto b = branchList.begin(); b != branchList.end(); ++b) {
    (*b)->set_new_substitutions();
  }

  traverse_find_ancestral_sequences();

  return(false);
}

bool Tree::step_through_MSAs() {
  // This is not going to work anymore as the tree resampling algorithm has changed.
  MSA->step_to_next_MSA();
  return(false);
}

// Record State data.
void Tree::initialize_output_streams() {
  files.add_file("tree", env.get("tree_out_file"), IOtype::OUTPUT);
  tree_out = files.get_ofstream("tree");
  files.add_file("substitutions", env.get("substitutions_out_file"), IOtype::OUTPUT);
  substitutions_out = files.get_ofstream("substitutions");

  substitutions_out << "Branch\tSubstitutions" << endl;
}

void Tree::record_tree() {
  /*
   * Records the tree topology with all node names and branch segments.
   */
  tree_out << root->toString();
}

void Tree::record_state(int gen, double l) {
  MSA->saveToFile(gen, l);
  //RecordSubtreeState();
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

void Tree::print_parameters() {
  SM->printParameters();	
}



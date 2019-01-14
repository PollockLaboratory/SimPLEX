#include "Tree.h"

#include <iostream>
#include <sstream> // For ostringstream
#include <cmath> // for floor and pow

#include "Sequence.h"
#include "Environment.h"
#include "IO.h"
#include "utils.h"

extern double Random();
extern Environment env;
extern IO::Files files;

using namespace std;

ofstream Tree::substitutions_out;
ofstream Tree::tree_out;

// Tree constructor.
Tree::Tree() {   
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

  // SM->clear_locations();
  for(auto b = branchList.begin(); b != branchList.end(); ++b) {
  	(*b)->update();
  }

  find_substitution_counts();
  InitializeOutputStreams();	

  std::cout << std::endl;
}

// Debug tools.
void Tree::printBranchList() {
  std::cout << "Printing Branch list. Size: " << branchList.size() << std::endl;
  for(std::list<BranchSegment*>::iterator it = branchList.begin(); it != branchList.end(); ++it) {
    BranchSegment* b = *it;
    std::cout << "Branch: " << b->distance;
    std::vector<substitution> subs = b->subs;
    for(auto s = subs.begin(); s != subs.end(); ++s) {
      std::cout << " " << MSA->decodeChar((*s).anc) << (*s).pos << MSA->decodeChar((*s).dec);
    }
    std::cout << std::endl;
  }
}

void Tree::printNodeList() {
  std::cout << "Printing Node list. Size:  " << nodeList.size() << std::endl;
  for(auto it = nodeList.begin(); it != nodeList.end(); ++it) {
    std::cout << "Node: " << (*it)->name << std::endl;
  }
}

void Tree::printParameters() {
  SM->printParameters();	
}

void Tree::printCounts() {
  std::cout << "Printing counts:"  << std::endl;
  std::cout << branchList.size() << std::endl;
  std::cout << nodeList.size() << std::endl;
  for(auto it = substitution_counts.begin(); it != substitution_counts.end(); ++it) {
    std::cout << "Distance: " << it->first << " 0: " << (it->second).first << " 1: " << (it->second).second << std::endl;
  }
}

// Sampling and likelihood.
void Tree::find_substitution_counts() {
  /*
   * Finds the key stats for the likelihood function.
   */

  SM->get_counts();

  substitution_counts = {};

  for(auto it = branchList.begin(); it != branchList.end(); ++it) {
    BranchSegment* b = *it;
    if(substitution_counts.find(b->distance) == substitution_counts.end()) {
      // Branch class does NOT already exists.
      substitution_counts[b->distance] = std::make_pair(b->num0subs, b->num1subs);
    } else {
      // Branch class does already exist.
      std::pair<int, int> c = substitution_counts[b->distance];
      substitution_counts[b->distance] = std::make_pair(b->num0subs + c.first, b->num1subs + c.second);
    }
  }
}

double Tree::calculate_likelihood() {
  double l_waiting = 0.0;
  float t;
  int num0subs;
  int num1subs;

  // Waiting times - this doesn't have to be calculated everytime.
  for(auto it = substitution_counts.begin(); it != substitution_counts.end(); ++it) {
    t = it->first;
    num0subs = it->second.first;
    num1subs = it->second.second;
    l_waiting += num0subs * log(1/(1 + u*t)) + num1subs * log(t/(1 + u*t));
  }

  double l_subs = SM->get_substitution_logLikelihood();

  // double l_old = 0.0;
  
  // Substitutions.
  // float r = 0;
  // for(auto b = branchList.begin(); b != branchList.end(); ++b) {
  // BranchSegment* branch = *b;
  // std::vector<substitution> subs = branch->subs;
  // for(int i = 0; i < subs.size(); i++) {
  //    if(subs[i].pos != -1) {
  //    r = branch->get_rate(i, subs[i].dec);
  //	if(isnan(log(r))) {
  //    std::cout << "Hit a nan rate: " << r << " " << log(r) << std::endl;
  //}
  //l_old += log(r);
  //   }
  // }
  //}

  return(l_waiting+l_subs);
}

void Tree::InitializeOutputStreams() {
  files.add_file("tree", env.get("tree_out_file"), IOtype::OUTPUT);
  tree_out = files.get_ofstream("tree");
  files.add_file("substitutions", env.get("substitutions_out_file"), IOtype::OUTPUT);
  substitutions_out = files.get_ofstream("substitutions");

  substitutions_out << "Branch\tSubstitutions" << endl;
}

//
// SAMPLING TREE PARAMETERS.
//

bool Tree::sample() {
  bool s = (this->*treeSamplingMethod)();

  // SM->clear_locations();
  // Update branch list - new substitutions.
  for(auto b = branchList.begin(); b != branchList.end(); ++b) {
    (*b)->update();
  }

  find_substitution_counts();

  return(s);
}

// Two options for changing ancestral states.
//

bool Tree::sample_ancestral_states() {
  // Find Random Node to start sampling.
  double r = Random();
  float s = 1.0 / nodeList.size();
  int i = 0;
  while(r > (i+1)*s) {
		i++;
  }

  // std::cout << "Starting node: " << i << " " <<  std::endl;
  nodeList[i]->sample();

  for(auto n = nodeList.begin(); n != nodeList.end(); ++n) {
    (*n)->sampled = false;
  }
  return(false);
}

bool Tree::step_through_MSAs() {
  MSA->step_to_next_MSA();
  return(false);
}

// Record State.
void Tree::RecordTree() {
  tree_out << root->toString();
}

void Tree::RecordState(int gen, double l) {
  MSA->saveToFile(gen, l);
  //RecordSubtreeState();
  //AddGenerationEndIndicatorsToOutputFiles();
}


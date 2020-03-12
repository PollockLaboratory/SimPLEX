#include "Tree.h"

#include <iostream>
#include <sstream> // For ostringstream
#include <cmath> // for floor and pow
#include <unordered_set>

#include "../Sequence.h"
#include "../../Environment.h"
#include "../../IO/Files.h"

extern double Random();
extern Environment env;
extern IO::Files files;

// Tree constructor.
Tree::Tree() {   
}

// Tree Initialize using seqs and states
void Tree::Initialize(IO::RawTreeNode* raw_tree, SequenceAlignment* &MSA, SubstitutionModel* &SM) {
  // Configuration
  max_seg_len = env.get<double>("TREE.max_segment_length");

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
  SM->organizeRateVectors(MSA->numCols(), SM->get_states()->n);

  // Proxys are created to correctly create root node.
  BranchSegment* proxyBranch = new BranchSegment(0.0);
  TreeNode* proxyNode = new TreeNode("Proxy");
  root = createTreeNode(raw_tree, proxyNode, proxyBranch);
  delete proxyBranch;
  delete proxyNode;

  configureSequences(root);

  for(auto t = nodeList.begin(); t != nodeList.end(); ++t) {
    if((*t)->isTip()) {
      tipList.push_back(*t);
    }
  }

  // Mark gaps.
  identify_gaps();
  
  // Initial sample to get counts.
  for(auto b = branchList.begin(); b != branchList.end(); ++b) {
  	(*b)->update();
  }

  //Setup output.
  files.add_file("tree_out", env.get<std::string>("OUTPUT.tree_out_file"), IOtype::OUTPUT);

  files.add_file("substitutions_out", env.get<std::string>("OUTPUT.substitutions_out_file"), IOtype::OUTPUT);
  files.write_to_file("substitutions_out", "I,GEN,LogL,Ancestral,Decendant,Substitutions\n");
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

SubstitutionModel* Tree::get_SM() {
  return(SM);
}

// Sampling and likelihood.
const std::list<BranchSegment*> Tree::get_branches() {
  return(branchList);
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

void Tree::identify_gaps() {
  for(auto n = nodeList.begin(); n != nodeList.end(); ++n) {
    (*n)->sampledp = false;
  }

  std::queue<TreeNode*> nodes = {};

  // Add tip nodes to starting nodes.
  for(auto t = tipList.begin(); t != tipList.end(); ++t) {
    (*t)->sampledp = false;
    nodes.push(*t);
  }

  while(not nodes.empty()) {
    if(not nodes.front()->sampledp) {
      nodes.push(nodes.front()->set_gaps());
    }

    nodes.pop();
  }

  // Reset all nodes such that the sampled flag is false.
  for(auto n = nodeList.begin(); n != nodeList.end(); ++n) {
    (*n)->sampledp = false;
  }
}

//
// SAMPLING TREE PARAMETERS.
//

sample_status Tree::sample(const std::list<int>& positions) {
  sample_status s = (this->*treeSamplingMethod)(positions);

  // Update branch list - new substitutions.
  for(auto b = branchList.begin(); b != branchList.end(); ++b) {
    (*b)->update();
  }

  return(s);
}

//
// Two options for changing ancestral states.
//

// Option 1.
sample_status Tree::sample_ancestral_states(const std::list<int>& positions) {
  /*
   * Both recalculates substitution events and recalculates the ancestral sequences.
   */
 
  // Reset all nodes such that the sampled flag is false.
  for(auto n = nodeList.begin(); n != nodeList.end(); ++n) {
    (*n)->sampledp = false;
  }

  std::queue<TreeNode*> nodes = {};
  
  // Add tip nodes to nodes.
  for(auto t = tipList.begin(); t != tipList.end(); ++t) {
    (*t)->sampledp = false;
    nodes.push(*t);
  }

  while(not nodes.empty()) {
    if(not nodes.front()->sampledp) {
      nodes.push(nodes.front()->calculate_state_probabilities(positions));
    }

    nodes.pop();
  }

  root->pick_sequences(positions);

  // Reset all nodes such that the sampled flag is false.
  for(auto n = nodeList.begin(); n != nodeList.end(); ++n) {
    (*n)->sampledp = false;
  }

  return(sample_status({false, true, true}));
}

// Option 2.
sample_status Tree::step_through_MSAs(const std::list<int>& positions) {
  // This is not going to work anymore as the tree resampling algorithm has changed.
  MSA->step_to_next_MSA();
  return(sample_status({false, true, true}));
}

// Record State data.
void Tree::record_tree() {
  /*
   * Records the tree topology with all node names and branch segments.
   */
  files.write_to_file("tree_out", root->toString());
}

void Tree::record_substitutions(int gen, double l) {
  static int index = -1;
  index++;

  std::ostringstream buffer;
  for(auto it = branchList.begin(); it != branchList.end(); ++it) {
    buffer << index << "," << gen << "," << l << ",";
    buffer << (*it)->ancestral->name << "," << (*it)->decendant->name << ",[ ";
    std::vector<Substitution> subs = (*it)->get_substitutions();
    for(unsigned int i = 0; i < subs.size(); i++) {
	TreeNode* anc = (*it)->ancestral;
	TreeNode* dec = (*it)->decendant;
	if(subs[i].occuredp == true) {
	  buffer << (*it)->ancestral->state_at_pos(i) << i << (*it)->decendant->state_at_pos(i) << " ";
	}
    }
    buffer << "]\n";
  }
  files.write_to_file("substitutions_out", buffer.str());
}

void Tree::record_state(int gen, double l) {
  MSA->saveToFile(gen, l);
  record_substitutions(gen, l);
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

// Ancestral States Parameter

AncestralStatesParameter::AncestralStatesParameter() : SampleableComponent("AncestralStates") {
  tree = new Tree();
  n_samples =  env.get<int>("MCMC.tree_sample_n_positions");
}

void AncestralStatesParameter::print() {
  std::cout << "AncestralStates" << std::endl;
}

std::string AncestralStatesParameter::get_type() {
  return("TREE_PARAMETER");
}

std::list<int> random_positions(int s_length, int n) {
  std::list<int> positions = {};
  if(n > s_length) {
    for(int i = 0; i < s_length; i++) {
      positions.push_back(i);
    }
  } else {
    // Very inefficient.
    int i;
    while(positions.size() < n) {
      i = rand() % s_length;
      positions.push_back(i);
      positions.unique();
    }
  }
  return(positions);
}

sample_status AncestralStatesParameter::sample() {
  // Pick positions.
  std::list<int> positions = random_positions(env.n, n_samples);
  return(tree->sample(positions));
}

void AncestralStatesParameter::undo() {
  std::cout << "Error: TreeParameter update cannot be undone." << std::endl;
  exit(EXIT_FAILURE);
}

void AncestralStatesParameter::fix() {
}

void AncestralStatesParameter::refresh() {
}

void AncestralStatesParameter::Initialize(IO::RawTreeNode* raw_tree, IO::RawMSA* &raw_msa, SubstitutionModel* &SM) {
  const States* states = SM->get_states();
  SequenceAlignment* MSA = new SequenceAlignment(states);
  MSA->Initialize(raw_msa);

  tree->Initialize(raw_tree, MSA, SM);
  tree->record_tree();
}

Tree* AncestralStatesParameter::get_tree_ptr() {
  return(tree);
}

double AncestralStatesParameter::record_state(int gen, double l) {
  tree->record_state(gen, l);
  return(0.0);
}

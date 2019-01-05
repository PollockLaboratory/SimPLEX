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
void Tree::connectNodes(TreeNode* &ancestralNode, BranchSegment* &ancestralBP, TreeNode* &decendantNode, float distance) {
	/* This function is responsible for connecting two Nodes with branch segments.
	 * This includes braking up long branches into smaller branch segment when needed and creating new internal branch nodes.
	 * Arguments:
	 * 	ancestralNode - pointer to the ancestral node.
	 * 	ancestralBP - pointer of the ancestral node that points to the decendant node. (This is needed to distinguish between the left and right branch pointers of ancestral node.)
	 * 	decendantNode - pointer to the decendantNode.
	 * 	distance - the branch length.
	 */
	std::pair<BranchSegment*, BranchSegment*> intermediateBranches = splitBranch(distance);
	BranchSegment* newBranchTop = intermediateBranches.first;
	BranchSegment* newBranchBottom = intermediateBranches.second;

	decendantNode->up = newBranchBottom;
	ancestralBP = newBranchTop;

	newBranchTop->ancestral = ancestralNode;
	newBranchBottom->decendant = decendantNode;
}

TreeNode* Tree::createTreeNode(IO::RawTreeNode* raw_tree, TreeNode* &ancestralNode, BranchSegment* &ancestralBP) {
	TreeNode* newTreeNode = new TreeNode(raw_tree);
	connectNodes(ancestralNode, ancestralBP, newTreeNode, raw_tree->distance);

	if(raw_tree->left != 0) { 
		createTreeNode(raw_tree->left, newTreeNode, newTreeNode->left);
	}

	if(raw_tree->right != 0) {
		createTreeNode(raw_tree->right, newTreeNode, newTreeNode->right); }

	return(newTreeNode);
}

// Configure sequences.

void Tree::configureSequences(TreeNode* n) {
	/*
	 * Traverses the tree attaching sequences to nodes.
	 * Also adds all the Nodes and Branch segments to their corresponding lists.
	 */

	nodeList.push_back(n);
	n->MSA = MSA;
	n->SM = SM;

	if(MSA->taxa_names_to_sequences.count(n->name)) {
		n->sequence = &(MSA->taxa_names_to_sequences.at(n->name));
	} else {
		if(n->isTip()){
			std::cerr << "Error: Missing sequence for \"" << n->name << "\"." << std::endl;
			exit(EXIT_FAILURE);
		}
		MSA->add(n->name);
		n->sequence = &(MSA->taxa_names_to_sequences.at(n->name));
	}

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

	// Fill missing sequences/
	if(n->isTip()) {
		// Node is tip.
	} else if(n->left != 0 and n->right == 0) {
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

void Tree::configureRateVectors() {
	BranchSegment* branch;
	std::vector<int> seq; 
	for(auto it = branchList.begin(); it != branchList.end(); ++it) {
		branch = (*it);
		seq = *(branch->ancestral->sequence);
		for(int i = 0; i < seq.size(); i++) {
			// Don't add rate vector if gap.
			if(seq[i] == -1) {
				branch->set_rate_vector(i);
			} else {
				RateVector* rv = SM->selectRateVector(seq[i]);
				branch->set_rate_vector(i, rv);	
			}
		}
	}
}

// Tree Initialize using seqs and states
void Tree::Initialize(IO::RawTreeNode* raw_tree, SequenceAlignment* &MSA, SubstitutionModel* &SM) {
	std::cout << "Creating MCMC tree structure." << std::endl;

	max_seg_len = env.get_float("max_segment_length");
	u = env.get_float("uniformization_constant");
	
	splitBranch = pickBranchSplitAlgorithm();
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
	std::cout << "Attaching rate vectors to tree." << std::endl;
	configureRateVectors();

	for(auto b = branchList.begin(); b != branchList.end(); ++b) {
		(*b)->updateStats();
	}

	findKeyStatistics();
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

// Sampling and likelihood.
void Tree::findKeyStatistics() {
	/*
	 * Finds the key stats for the likelihood function.
	 */

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
	double l = 0.0;

	// Waiting times.
	for(auto it = substitution_counts.begin(); it != substitution_counts.end(); ++it) {
		float t = it->first;
		int num0subs = it->second.first;
		int num1subs = it->second.second;
		l += num0subs * log(1/(1 + u*t)) + num1subs * log(t/(1 + u*t));
	}

	// Substitutions.
	float r = 0;
	for(auto b = branchList.begin(); b != branchList.end(); ++b) {
		BranchSegment* branch = *b;
		std::vector<substitution> subs = branch->subs;
		for(int i = 0; i < subs.size(); i++) {
			if(subs[i].pos != -1) {
				r = branch->get_rate(i, subs[i].dec);
				if(isnan(log(r))) {
					std::cout << "Hit a nan rate: " << r << " " << log(r) << std::endl;
				}
				l += log(r);
			}
		}
	}
	return(l);
}

double Tree::partial_calculate_likelihood() {
	//std::cout << "Partial Likelihood Calculation - tree" << std::endl;
	double deltaLogL = 0.0;
	std::list<AbstractValue*> l = SM->get_current_parameters();
	for(auto it = l.begin(); it != l.end(); ++it) {
		int dec = (*it)->state;
		std::set<bpos> locs = (*it)->rv->get_locations();
		for(auto jt = locs.begin(); jt != locs.end(); ++jt) {
			BranchSegment* b = (*jt).branch;
			substitution s = b->subs[(*jt).pos];
			if(dec == s.dec) {
				deltaLogL += log((*it)->getOldValue()/(*it)->getValue());	
			}
		}
	}

	return(deltaLogL);
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

void SampleNode(TreeNode* n) {
	if(n->sampled == false) {
		n->sampled = true;
		n->sampleSequence();

		if(n->left != 0) {
			SampleNode(n->left->decendant);
		}

		if(n->right != 0) {
			SampleNode(n->right->decendant);
		}
		
		if(n->up != 0) {
			SampleNode(n->up->ancestral);
		}
	}
}

bool Tree::SampleParameters() {
	double r = Random();
	float s = 1.0 / nodeList.size();
	int i = 0;
	while(r > (i+1)*s) {
		i++;
	}

	// std::cout << "Starting node: " << i << " " <<  std::endl;
	SampleNode(nodeList[i]);

	for(auto n = nodeList.begin(); n != nodeList.end(); ++n) {
		(*n)->sampled = false;
	}
	
	// Update branch list - new substitutions.
	for(auto b = branchList.begin(); b != branchList.end(); ++b) {
		(*b)->updateStats();
	}

	findKeyStatistics();

	return(false);
}

// Record State.
void Tree::RecordTree() {
  std::cout << "Recording Tree." << std::endl;
  std::cout << root->name << std::endl;
  tree_out << root->toString();
}

void Tree::RecordState(int gen, double l) {
  MSA->saveToFile(gen, l);
  //RecordSubtreeState();
  //AddGenerationEndIndicatorsToOutputFiles();
}


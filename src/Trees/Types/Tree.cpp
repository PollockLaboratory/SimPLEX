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

int Tree::num_trees = 0;

ofstream Tree::substitutions_out;
ofstream Tree::tree_out;

// Tree constructor.
Tree::Tree() {   
	std::cout << "Creating Basic tree." << std::endl;
	is_constant = true;
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
		b->rates.resize(seqLen);
		branchList.push_back(b);
		configureSequences(b->decendant);
	}

	if(n->right != 0) {
		BranchSegment* b = n->right;
		b->rates.resize(seqLen);
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
	std::cout << "Configuring rate vectors." << std::endl;
	BranchSegment* branch;
	std::vector<int> s; 
	for(auto it = branchList.begin(); it != branchList.end(); ++it) {
		branch = (*it);
		s = *(branch->ancestral->sequence);
		for(int i = 0; i < s.size(); i++) {
			// Don't add rate vector if gap.
			if(s[i] == -1) {
				branch->rates[i] = NULL;
			} else {
				branch->rates[i] = SM->selectRateVector(s[i]);	
			}
		}
	}
}

// Tree Initialize using seqs and states
void Tree::Initialize(IO::RawTreeNode* raw_tree, SequenceAlignment* &MSA, SubstitutionModel* &SM) {
	std::cout << "INITIALIZING BASIC TREE." << std::endl;
	max_seg_len = env.get_float("max_segment_length");
	std::cout << "Max segment length: " << max_seg_len << std::endl;

	u = env.get_float("uniformization_constant");

	splitBranch = pickBranchSplitAlgorithm();
	this->MSA = MSA;
	seqLen = MSA->numCols();
	this->SM = SM;

	// Proxys are created to correctly create root node.
	BranchSegment* proxyBranch = new BranchSegment(0.0);
	TreeNode* proxyNode = new TreeNode("Proxy");
	TreeNode* root = createTreeNode(raw_tree, proxyNode, proxyBranch);
	delete proxyBranch;
	delete proxyNode;

	configureSequences(root);
	configureRateVectors();

	for(auto b = branchList.begin(); b != branchList.end(); ++b) {
		(*b)->updateStats();
	}

	findKeyStatistics();
	InitializeOutputStreams();
	
	// printNodeList();
	// printBranchList();
}

// Debug tools.

void Tree::printBranchList() {
	std::cout << "Printing Branch list. Size: " << branchList.size() << std::endl;
	for(std::list<BranchSegment*>::iterator it = branchList.begin(); it != branchList.end(); ++it) {
		BranchSegment* b = *it;
		std::cout << "Branch: " << b->distance;
		auto subs = b->subs;
		for(std::list<substitution>::iterator s = subs.begin(); s != subs.end(); ++s) {
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
	for(auto b = branchList.begin(); b != branchList.end(); ++b) {
		BranchSegment* branch = *b;
		for(auto s = branch->subs.begin(); s != branch->subs.end(); ++s) {
			substitution sub = *s;
			l += log(branch->rates[sub.pos]->rates[sub.dec]->getValue());
		}
	}

	return(l);
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
		// std::cout << "Sampling tree node: " << n->name << std::endl;
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
void Tree::RecordState(int gen, double l) {
	MSA->saveToFile(gen, l);
	//RecordSubtreeState();
	//AddGenerationEndIndicatorsToOutputFiles();
}


#include "Tree.h"

#include <iostream>
#include <sstream> // For ostringstream
#include <cmath> // for floor and pow
#include "Environment.h"
#include "IO.h"
#include "utils.h"

extern double Random();
extern Environment env;
extern IO::Files files;

using namespace std;

int Tree::num_trees = 0;

ofstream Tree::substitutions_out;
ofstream Tree::sequences_out;
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

	std::cout << "Configure: " << n->name << std::endl;
	nodeList.push_back(n);
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
		std::cout << "Tip Here: " << n->name << std::endl;
	} else if(n->left != 0 and n->right == 0) {
		// Internal Continous.
		std::cout << "Internal Continous Here: " << n->name << std::endl;
		TreeNode* dsNode = n->left->decendant; // ds = downstream.
		*(n->sequence) = *(dsNode->sequence);
	} else {
		// Root or internal branch.
		std::cout << "Internal Splitting Here: " << n->name << std::endl;
		vector<int> dsNodeLseq = *(n->left->decendant->sequence);
		vector<int> dsNodeRseq = *(n->right->decendant->sequence);
		vector<int> p = MSA->findParsimony(dsNodeLseq, dsNodeRseq);
		*(n->sequence) = p;

		// Add substituions to branches.
		n->left->subs = MSA->findSubstitutions(p, dsNodeLseq);
		n->right->subs = MSA->findSubstitutions(p, dsNodeRseq);
	}
}

// Tree Initialize using seqs and states
void Tree::Initialize(IO::RawTreeNode* raw_tree, SequenceAlignment* &MSA) {
	std::cout << "INITIALIZING BASIC TREE." << std::endl;
	max_seg_len = env.get_float("max_segment_length");
	std::cout << "Max segment length: " << max_seg_len << std::endl;

	splitBranch = pickBranchSplitAlgorithm();
	this->MSA = MSA;

	// Proxys are created to correctly create root node.
	BranchSegment* proxyBranch = new BranchSegment(0.0);
	TreeNode* proxyNode = new TreeNode("Proxy");
	TreeNode* root = createTreeNode(raw_tree, proxyNode, proxyBranch);
	delete proxyBranch;
	delete proxyNode;

	configureSequences(root);
	printNodeList();
	printBranchList();

	MSA->print();
	//InitializeSequences(taxa_names_to_sequences);
	InitializeOutputStreams();
	//RecordState();
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
	for(std::list<TreeNode*>::iterator it = nodeList.begin(); it != nodeList.end(); ++it) {
		std::cout << "Node: " << (*it)->name << std::endl;
	}

}

// Sampling and likelihood.

double Tree::calculate_likelihood() {
	return(0.0);
}

// bool Tree::IsRoot() const { return up == NULL;  }
// bool Tree::IsSubtree() const { return (left != NULL and right != NULL); }
// bool Tree::IsLeaf() const { return (left == NULL and right == NULL); }

//STP: I have chosen the name "segment" to describe a node that exists only
// to subdivide a branch.
// Segments will have a node to the left and not to the right
// There are no segments yet.
//bool Tree::IsSegment() const {
//	return (left != NULL and right == NULL);
//}

/**
 *
 * Now one can initialize internal node sequences
 *
 */

// void Tree::InitializeSequences(
// 		map<string, vector<int> > taxa_names_to_sequences) {
// 	if (IsSubtree()) {
// 		left->InitializeSequences(taxa_names_to_sequences);
// 		right->InitializeSequences(taxa_names_to_sequences);
// 
// 		sequence.resize(left->sequence.size());
// 		SampleSequence();
// 
// 	} else if (IsLeaf()) {
// 		if (taxa_names_to_sequences.find(name)
// 				== taxa_names_to_sequences.end()) {
// 			cerr << "Could not find a sequence for " << name << endl;
// 			exit(-1);
// 		}
//
//		sequence = taxa_names_to_sequences[name];
//	}
//}

void Tree::InitializeOutputStreams() {
	files.add_file("tree", env.get("tree_out_file"), IOtype::OUTPUT);
	tree_out = files.get_ofstream("tree");
	files.add_file("sequences", env.get("sequences_out_file"), IOtype::OUTPUT);
	sequences_out = files.get_ofstream("sequences");
	files.add_file("substitutions", env.get("substitutions_out_file"), IOtype::OUTPUT);
	substitutions_out = files.get_ofstream("sequences");

	substitutions_out << "Branch\tSubstitutions" << endl;
}

//
// SAMPLING TREE PARAMETERS.
//

void Tree::SampleParameters() {
	//std::cout << "Sampling tree parameters" << std::endl;
	//SampleSubtreeParameters();
}

/**
 * This must be post-order in order to update the ancestors from the
 * descendants. Need to sample the descendants' sequences before the
 * ancestral sequences can be sampled.
 *
 */
//void Tree::SampleSubtreeParameters() {
//	if (IsSubtree()) {
//		left->SampleSubtreeParameters();
//		right->SampleSubtreeParameters();
//	}
//	SampleSequence();
//	SampleDistance();
//}

/**
 * There are many different ancestral state sampling algorithms.
 *
 */

//void Tree::SampleSequence() {
//	if (not IsLeaf()) {  MetropolisHastingsStateSampling(); }
//}

/**
 *
 * This is the simplest ancestral state sampling method I could think of.
 * Next, one might imagine allowing the ancestor to
 * be the same state as the descendant but choose the left or the right with
 * the probability of the substitution required. For example, then low
 * probability substitutions would not be required and the MCMC would sample
 * the high probability substitutions more often than the low probability
 * substitutions.
 * p(left) = p(substitution required if ancestor state is left state) / normalize
 * p(right) = p(substitution required if ancestor state is right state) / normalize
 *
 * That would be Gibbs sampling because it samples "from the posterior" meaning
 * it takes into account the relative likelihood values of all the states
 * to inform sampling.
 *
 * This sampling method has a strong prior on all nodes and all sites. The
 * prior essentially is produced by the states at the leaves because this site
 * at this node can only take the values of a state at a leaf.
 *
 * This does follow detailed balance but has an overpowering prior and so will
 * produce results biased using that prior.
 *
 *
 */

//void Tree::DescendentStateSampling() {
//	for (int site = 0; site < sequence.size(); site++) {
//		if (Random() < 0.5) {
//			sequence.at(site) = left->sequence.at(site);
//		} else {
//			sequence.at(site) = right->sequence.at(site);
//		}
//	}
//}

/*  uniform prior for sampling thestates at all sites for a node.
 The sample is not influenced by data. This method should mix poorly.
 */
//void Tree::MetropolisHastingsStateSampling() { // this is not MH
//	for (int site = 0; site < sequence.size(); site++) {
//		int state_integer = std::floor(Random() * states.size());
//		sequence.at(site) = state_integer;
//	}
//}

/**
 * This is a simple sampling method to take the previous value into account.
 * The new value is allowed to be 50% greater or lower than the previous value.
 * This allows the new value to sample larger and smaller values.
 *
 *
 *
 * This method, however, does not follow detailed balance because the
 * transition rate (and probability in this case) forward and back are not the
 * same. If one jumps down by 50%, one cannot jump back up to the original
 * value because it is 200% of the current value.
 * Taking a constant  sized step in either direction would follow detailed
 * balance but then I would have to choose a step size and I don't want to do
 * that.
 *
 * Or I could choose a small step size like up by 10/9 or down by 9/10 and that
 * would be detail balanced because the probability of going there and back
 * are equal.
 *
 * I want the value to be able to wander up or down equally so I determined
 * the map from a random value between 0 and 1 to correspond to a decrease to
 * 66% or an increase to 150% respectively.
 *
 * x = 0, y = 2/3
 * x = 1/2, y = 1
 * x = 1, y = 1.5
 *
 * y = 4/3 x ^ 2 + 2/3
 *
 * This method follows detailed balance because the model has equal transition
 * probabilities from one value to a different value and back again.
 *
 */

// simpler to make step size, choose random
// could also sometimes choose smaller or larger step size
//void Tree::SampleDistance() {
//	if (not is_constant) {
//		//distance = distance * (4.0 / 3.0 * std::pow(Random(), 2) + 2.0 / 3.0);
//		if (Random() < 0.5) { distance *= 9.0 / 10.0;
//		} else { distance *= 10.0 / 9.0; }
//	}
//}

void Tree::RecordState() {
	//RecordSubtreeState();
	//AddGenerationEndIndicatorsToOutputFiles();
}

//STP: RecordState should be recursive so long as you traverse the tree post-
// order. Newick tree format is a post-order format. Order does not matter for
// the sequences and substitutions.
//void Tree::RecordSubtreeState() {
//	if (IsSubtree()) {
//		tree_out << "(";
//		left->RecordSubtreeState();
//		tree_out << ",";
//		right->RecordSubtreeState();
//		tree_out << ")";
//	}
//	RecordToTreefile();
//	RecordSequence();
//	RecordSubstitutions();
//}

/* Need to remember to print tree at least once to attach the names of the internal nodes to the tree. This needs a better name.  */
//void Tree::RecordToTreefile() {
//	tree_out << name << ":" << distance;
//}

//STP: Recording the sequences requires decoding the sequences from integers to
// amino acids or nucleotides (chars).
//void Tree::RecordSequence() {
//	sequences_out << ">" << name << endl;
	
//	for (int site = 0; site < sequence.size(); site++) {
//		sequences_out << states[sequence.at(site)];
//	}
//	sequences_out << endl;
//}

//STP: When are the substitutions sampled? Only when we sample the internal
// sequences. Then we could determine the substitutions once then. I don't
// want to figure out the substitutions again here.

//void Tree::RecordSubstitutions() {
	// At the moment, the substitutions are calculated at record time. (See
	// function "_printsitespecificlist". Does this make the most sense?
	// Can't we calculate the substitutions when we resample the ancestral
	// states? We might sample the ancestral states much more often than
	// recording the substitutions. For now I will calculate the substitutions
	// here.

	/*
	 * It would be better if this function could record the substitutions from
	 * its parent to itself but that would require knowing who its parent is.
	 * Instead, I will make this function print the substitutions from this
	 * parent to its children. That is, the substitutions attached to the
	 * children nodes.
	 * This is better because the parent knows who its children are while the
	 * children don't know who their parents are.
	 *
	 * This also avoids the problem of calling this function on the root node.
	 * The root node shouldn't know it is the root node.
	 */

//	if (IsSubtree()) {
//		RecordChildSubstitutions(left);
//		RecordChildSubstitutions(right);
//	}

//}

/**
 * The name of this function does not make it obvious what the argument is.
 *
 * Since I do this for both children and I didn't want to repeat myself
 * (DRY principle)
 * This function requires the global Data data for the reverse state index.
 * It might make reverseStateIndex a global since the tree needs to use it.
 * But only the tree needs to use it; no other part needs the information.
 * Should the tree be able to decode itself?
 *
 */
//void Tree::RecordChildSubstitutions(Tree* child) {
//	substitutions_out << child->name;
//	for (int site = 0; site < sequence.size(); site++) {
//		if (sequence.at(site) != child->sequence.at(site)) {
//			//Note: This prints out one-indexed sequence positions.
//			substitutions_out << "\t" << states[sequence.at(site)] << (site + 1)
//					<< states[child->sequence.at(site)];
//		}
//	}
//	substitutions_out << endl;
//}

// simplify explanation here
/* Add final semicolon only at root node. But root doesn't know it is root, so call this function after the root has recorded itself.   */
//void Tree::AddGenerationEndIndicatorsToOutputFiles() {
//	tree_out << ";" << endl;
//	sequences_out << endl;
	// This is required because the root node is not a child of anything and so
	// its substitutions are never printed. This would not be here if a tree
	// could record its own substitutions.
//	substitutions_out << name << endl;
//	substitutions_out << endl;
//}

#include "Tree_B1.h"

#include <iostream>
#include <sstream> // For ostringstream
#include <cmath> // for floor and pow
#include "Options.h"

extern double Random();
extern Options options;
using namespace std;

// These are instantiated in Tree.cpp
//int Tree::number_of_trees = 0;

// See "Print yourself.odt" for rationale for why these are class static.
//ofstream Tree_B1::substitutions_out;
//ofstream Tree_B1::sequences_out;
//ofstream Tree_B1::tree_out;

Tree_B1::Tree_B1() {
	// Is this considered doing work? Is this exception safe?
	/* All this is done in the tree constructor
	
	id = number_of_trees;
	number_of_trees++;

	is_constant = false;

	distance = 0;

	left = NULL;
	right = NULL;
	up = NULL;
	*/
	/* if (options.debug)
		std::cout << "Created tree b1 with id " << id << std::endl; */
}

Tree_B1::Tree_B1(const Tree_B1& tree) {
	id = num_trees;
	num_trees++;

	is_constant = tree.is_constant;

	name = tree.name;
	distance = tree.distance;
	sequence = tree.sequence;

	states = tree.states;

	if (tree.IsSubtree()) {
		left = tree.left->Clone();
		left->up = this;
		right = tree.right->Clone();
		right->up = this;
	} else if (tree.IsSegment()){
		left = tree.left->Clone();
		left->up = this;
		right = NULL;
	} else {
		left = NULL;
		right = NULL;
	}
}

//Copy-swap
Tree_B1& Tree_B1::operator=(Tree_B1 tree) {

	std::swap(id, tree.id);
	std::swap(is_constant, tree.is_constant);
	std::swap(name, tree.name);
	std::swap(distance, tree.distance);
	std::swap(sequence, tree.sequence);
	std::swap(left, tree.left);
	std::swap(right, tree.right);
	std::swap(up, tree.up);
	std::swap(states, tree.states);

	return *this;
}

// The Tree destructor is sufficient
/* Tree_B1::~Tree_B1() {
	
	if (IsSubtree()) {
		delete left;
		//left = NULL;
		delete right;
		// right = NULL;
	} else if (IsSegment()) {
		delete left;
		// left = NULL;
	} 
} */

Tree_B1* Tree_B1::Clone() {
//	std::cout << "cloning tree" << std::endl;
	return new Tree_B1(*this);
}

// This is the same as in Tree and might not need to be overridden
void Tree_B1::Initialize(map<string, vector<int> > taxa_names_to_sequences,
		vector<string> states) {
	this->states = states;
	ReadFromTreeFile();
	InitializeSequences(taxa_names_to_sequences);
	InitializeOutputStreams();
	RecordState();
}

void Tree_B1::ReadFromString(string tree_string) {
	// STP: Why is the constancy of the tree determined here?
	is_constant = options.constant_tree;
	
	ExtractDistance(tree_string);
	ExtractName(tree_string);
	
	if (distance > options.max_segment_length) {
		SegmentBranch();
	}

	if (tree_string.find('(') != string::npos
			and tree_string.find(')') != string::npos) {
		ReadSubtree(tree_string);
	}
}

void Tree_B1::SegmentBranch(){
	// At the end of this method, "this" must be the bottom segment and the top
	// segment's "up" must point to the current "this->up"
	int number_of_segments = std::ceil(distance / options.max_segment_length);
	float segment_length = distance / number_of_segments;
	
	std::cout << "This branch needs to be broken up " << name << std::endl
			  << "because the length is " << distance << std::endl;

	Tree* current_up = this->up;
	for (int i = 1; i <= number_of_segments; i++) {
		Tree* segment = new Tree_B1();
		segment->up = current_up;
		segment->distance = segment_length;
		segment->states = this->states;
		
		if (not IsRoot()){
			if (current_up == this->up and this->up->right == this) {
				current_up->right = segment;
			} else {
				current_up->left = segment;
			}
		}
		current_up = segment;
	}
	this->up = current_up;
	current_up->left = this;
	distance = segment_length;
}

void Tree_B1::ReadSubtree(string tree_string) {

	int tree_depth = 0;
	for (int position = 0; position < tree_string.size(); position++) {
		char character = tree_string.at(position);
		if (character == '(')
			tree_depth++;
		if (character == ')')
			tree_depth--;

		if (character == ',' and tree_depth == 1) {
			// The only difference is you make a new Tree_B1 rather than a new
			// Tree
			left = new Tree_B1();
			left->up = this;
			left->states = this->states;
			left->ReadFromString(tree_string.substr(1, position - 1));

			right = new Tree_B1();
			right->up = this;
			right->states = this->states;
			right->ReadFromString(
					tree_string.substr(position + 1,
							tree_string.size() - position - 2));
		}
	}
}


bool Tree_B1::IsSegment() const {
	return (left != NULL and right == NULL);
}

/**
 *
 * Now one can initialize internal node sequences
 *
 */
void Tree_B1::InitializeSequences(
		map<string, vector<int> > taxa_names_to_sequences) {
	
	if (IsSubtree()) {
		left->InitializeSequences(taxa_names_to_sequences);
		right->InitializeSequences(taxa_names_to_sequences);

		sequence.resize(left->sequence.size());
		SampleSequence();
	} else if (IsSegment()) {
		left->InitializeSequences(taxa_names_to_sequences);
		sequence.resize(left->sequence.size());
		SampleSequence();
	} else if (IsLeaf()) {
		if (taxa_names_to_sequences.find(name)
				== taxa_names_to_sequences.end()) {
			cerr << "Could not find a sequence for B1 " << name << endl;
			exit(-1);
		}

		sequence = taxa_names_to_sequences[name];
	}
}


void Tree_B1::SampleParameters() {
//	std::cout << "Sampling tree parameters" << std::endl;
	SampleSubtreeParameters();
}

/**
 * This must be post-order in order to update the ancestors from the
 * descendants. Need to sample the descendants' sequences before the
 * ancestral sequences can be sampled.
 *
 */
void Tree_B1::SampleSubtreeParameters() {
	if (IsSubtree()) {
		left->SampleSubtreeParameters();
		right->SampleSubtreeParameters();
	} else if (IsSegment()) {
		left->SampleSubtreeParameters();
	}
	SampleSequence();
	SampleDistance();
}



//STP: RecordState should be recursive so long as you traverse the tree post-
// order. Newick tree format is a post-order format. Order does not matter for
// the sequences and substitutions.
void Tree_B1::RecordSubtreeState() {
	if (IsSubtree()) {
		tree_out << "(";
		left->RecordSubtreeState();
		tree_out << ",";
		right->RecordSubtreeState();
		tree_out << ")";
	}
	if (IsSegment()) {
		tree_out << "(";
		left->RecordSubtreeState();
		tree_out << ")";
	}

	RecordToTreefile();
	RecordSequence();
	RecordSubstitutions();
}

//STP: When are the substitutions sampled? Only when we sample the internal
// sequences. Then we could determine the substitutions once then. I don't
// want to figure out the substitutions again here.
void Tree_B1::RecordSubstitutions() {
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

	if (IsSubtree()) {
		RecordChildSubstitutions(left);
		RecordChildSubstitutions(right);
	} else if (IsSegment()) {
		RecordChildSubstitutions(left);
	}

}

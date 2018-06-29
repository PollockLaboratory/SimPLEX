#include "Tree.h"

#include <iostream>
#include <sstream> // For ostringstream
#include <cmath> // for floor and pow
#include "../Options.h"
#include "../utils.h"

extern double Random();
extern Options options;
using namespace std;

int Tree::num_trees = 0;

// See "Print yourself.odt" for rationale for why these are class static.
ofstream Tree::substitutions_out;
ofstream Tree::sequences_out;
ofstream Tree::tree_out;

// Tree constructor.
Tree::Tree() {   
	std::cout << "Creating default tree." << std::endl;
	id = num_trees;  
	num_trees++;
	name = "Node_" + IdToString();
	is_constant = false;
	distance = 0;
	left = NULL; right = NULL; up = NULL;
}

// copy
Tree::Tree(const Tree& tree) { id = num_trees; num_trees++;
    name = tree.name;
    is_constant = tree.is_constant;
    distance = tree.distance;

    sequence = tree.sequence; states = tree.states; // not in constructor
	if (tree.IsSubtree()) {
		left = tree.left->Clone();  left->up = this;
		right = tree.right->Clone();  right->up = this;
	} else { left = NULL;  right = NULL;  }
}

// Tree destruction
Tree::~Tree() {
	if (left != NULL)  delete left;
	if (right != NULL) delete right;
}

// Tree cloning
Tree* Tree::Clone() {
	// std::cout << "cloning tree" << std::endl;
	return new Tree(*this);
}

// Tree Initialize using seqs and states
void Tree::Initialize(map<string, vector<int> > taxa_names_to_sequences, vector<string> states) {
	std::cout << "INITIALIZING TREE." << std::endl;
	this->states = states;
	ReadFromTreeFile();
	InitializeSequences(taxa_names_to_sequences);
	InitializeOutputStreams();
	RecordState();
}

void Tree::ReadFromTreeFile() {
	string treefile = options.get("tree_file");
	ifstream tree_in(treefile.c_str());

	if (not tree_in.good()) {
		cout << "Could not read tree file " << treefile << endl;
		exit(-1);
	}

	string tree_string;
	getline(tree_in, tree_string);

	//STP: Get rid of final semicolon and any spaces
	tree_string.erase(std::remove(tree_string.begin(), tree_string.end(), ';'),
			tree_string.end());
	tree_string.erase(std::remove(tree_string.begin(), tree_string.end(), ' '),
			tree_string.end());

	ReadFromString(tree_string);
}

void Tree::ReadFromString(string tree_string) {

	// STP: Why is the constancy of the tree determined here?
	is_constant = options.get_int("constant_tree");
	
	ExtractDistance(tree_string);
	ExtractName(tree_string);
	
	if (tree_string.find('(') != string::npos
			and tree_string.find(')') != string::npos) {
		ReadSubtree(tree_string);
	}
}

void Tree::ReadSubtree(string tree_string) {

	int tree_depth = 0;
	for (int position = 0; position < tree_string.size(); position++) {
		char character = tree_string.at(position);
		if (character == '(')
			tree_depth++;
		if (character == ')')
			tree_depth--;

		if (character == ',' and tree_depth == 1) {
			left = new Tree();
			left->states = this->states;
			left->up = this;
			left->ReadFromString(tree_string.substr(1, position - 1));

			right = new Tree();
			right->states = this->states;
			right->up = this;
			right->ReadFromString(tree_string.substr(position + 1,
							      tree_string.size() - position - 2));
		}
	}
}

std::string Tree::IdToString() {
	ostringstream id_str;
	id_str << id;
	return id_str.str();
}

void Tree::ExtractDistance(string& tree_string) {
	int last_colon_position = tree_string.find_last_of(':');

	if (last_colon_position == string::npos) {
		cerr << "Node " << name << " with id " << id
				<< " doesn't have a distance" << endl;
		exit(-1);
	}

	distance = atof(tree_string.substr(last_colon_position + 1).c_str());
	tree_string.resize(last_colon_position);
}

void Tree::ExtractName(string& tree_string) {
	int last_parens_position = tree_string.find_last_of(')');

	if (last_parens_position == string::npos) {
		name = tree_string;
	} else {
		name = tree_string.substr(last_parens_position + 1);
		tree_string.resize(last_parens_position + 1);
	}
}

bool Tree::IsRoot() const { return up == NULL;  }
bool Tree::IsSubtree() const { return (left != NULL and right != NULL); }
bool Tree::IsLeaf() const { return (left == NULL and right == NULL); }

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

void Tree::InitializeSequences(
		map<string, vector<int> > taxa_names_to_sequences) {
	if (IsSubtree()) {
		left->InitializeSequences(taxa_names_to_sequences);
		right->InitializeSequences(taxa_names_to_sequences);

		sequence.resize(left->sequence.size());
		SampleSequence();

	} else if (IsLeaf()) {
		if (taxa_names_to_sequences.find(name)
				== taxa_names_to_sequences.end()) {
			cerr << "Could not find a sequence for " << name << endl;
			exit(-1);
		}

		sequence = taxa_names_to_sequences[name];
	}
}

void Tree::InitializeOutputStreams() {
	tree_out.open(options.treeout.c_str());
	sequences_out.open(options.seqsout.c_str());
	substitutions_out.open(options.subsout.c_str());

	substitutions_out << "Branch\tSubstitutions" << endl;
}

// sampling tree parameters
void Tree::SampleParameters() { //	std::cout << "Sampling tree parameters" << std::endl;
	SampleSubtreeParameters();
}

/**
 * This must be post-order in order to update the ancestors from the
 * descendants. Need to sample the descendants' sequences before the
 * ancestral sequences can be sampled.
 *
 */
void Tree::SampleSubtreeParameters() {
	if (IsSubtree()) {
		left->SampleSubtreeParameters();
		right->SampleSubtreeParameters();
	}
	SampleSequence();
	SampleDistance();
}

/**
 * There are many different ancestral state sampling algorithms.
 *
 */

void Tree::SampleSequence() {
	if (not IsLeaf()) {  MetropolisHastingsStateSampling(); }
}

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

void Tree::DescendentStateSampling() {
	for (int site = 0; site < sequence.size(); site++) {
		if (Random() < 0.5) {
			sequence.at(site) = left->sequence.at(site);
		} else {
			sequence.at(site) = right->sequence.at(site);
		}
	}
}

/*  uniform prior for sampling thestates at all sites for a node.
 The sample is not influenced by data. This method should mix poorly.
 */
void Tree::MetropolisHastingsStateSampling() { // this is not MH
	for (int site = 0; site < sequence.size(); site++) {
		int state_integer = std::floor(Random() * states.size());
		sequence.at(site) = state_integer;
	}
}

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
void Tree::SampleDistance() {
	if (not is_constant) {
		//distance = distance * (4.0 / 3.0 * std::pow(Random(), 2) + 2.0 / 3.0);
		if (Random() < 0.5) { distance *= 9.0 / 10.0;
		} else { distance *= 10.0 / 9.0; }
	}
}

void Tree::RecordState() {
	RecordSubtreeState();
	AddGenerationEndIndicatorsToOutputFiles();
}

//STP: RecordState should be recursive so long as you traverse the tree post-
// order. Newick tree format is a post-order format. Order does not matter for
// the sequences and substitutions.
void Tree::RecordSubtreeState() {
	if (IsSubtree()) {
		tree_out << "(";
		left->RecordSubtreeState();
		tree_out << ",";
		right->RecordSubtreeState();
		tree_out << ")";
	}
	RecordToTreefile();
	RecordSequence();
	RecordSubstitutions();
}

/* Need to remember to print tree at least once to attach the names of the internal nodes to the tree. This needs a better name.  */
void Tree::RecordToTreefile() {
	tree_out << name << ":" << distance;
}

//STP: Recording the sequences requires decoding the sequences from integers to
// amino acids or nucleotides (chars).
void Tree::RecordSequence() {
	sequences_out << ">" << name << endl;
	
	for (int site = 0; site < sequence.size(); site++) {
		sequences_out << states[sequence.at(site)];
	}
	sequences_out << endl;
}

//STP: When are the substitutions sampled? Only when we sample the internal
// sequences. Then we could determine the substitutions once then. I don't
// want to figure out the substitutions again here.
void Tree::RecordSubstitutions() {
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
	}

}

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
void Tree::RecordChildSubstitutions(Tree* child) {
	substitutions_out << child->name;
	for (int site = 0; site < sequence.size(); site++) {
		if (sequence.at(site) != child->sequence.at(site)) {
			//Note: This prints out one-indexed sequence positions.
			substitutions_out << "\t" << states[sequence.at(site)] << (site + 1)
					<< states[child->sequence.at(site)];
		}
	}
	substitutions_out << endl;
}

// simplify explanation here
/* Add final semicolon only at root node. But root doesn't know it is root, so call this function after the root has recorded itself.   */
void Tree::AddGenerationEndIndicatorsToOutputFiles() {
	tree_out << ";" << endl;
	sequences_out << endl;
	// This is required because the root node is not a child of anything and so
	// its substitutions are never printed. This would not be here if a tree
	// could record its own substitutions.
	substitutions_out << name << endl;
	substitutions_out << endl;
}


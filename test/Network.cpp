#include "Network.h"

#include <sstream> // For ostringstream
#include "Options.h"
#include "Data.h"

extern double Random();

extern Options options;

using namespace std;

int Network::number_of_networks = 0;

// See "Print yourself.odt" for rational for why these are class static.
ofstream Network::substitutions_out;
ofstream Network::sequences_out;
ofstream Network::Network_out;

Network::Network() {
	id = number_of_networks;
	number_of_networks++;

	is_constant = false;

//	std::cout << "Network " << id << " constructed" << std::endl;

}

Network::Network(const Network& Network) {
	id = number_of_networks;
	number_of_networks++;

	is_constant = Network.is_constant;

	distances = Network.distances;

	integer_to_state = Network.integer_to_state;
}

void Network::Initialize(map<string, vector<int> > taxa_names_to_sequences,
		map<int, string> integer_to_state) {
	this->integer_to_state = integer_to_state;
	InitializeSequences(taxa_names_to_sequences);
	InitializeOutputStreams();
	RecordState();
	exit(1);
}

std::string Network::IdToString() {
	ostringstream id_str;
	id_str << id;
	return id_str.str();
}

void Network::InitializeSequences(map<string, vector<int> > taxa_names_to_sequences) {
	for (int i = 0; i < nodes.size(); i++) {
		nodes.at(i).sequence = taxa_names_to_sequences[nodes.at(i).name];
	}
}

void Network::InitializeOutputStreams() {
	Network_out.open(options.tree_out_file.c_str());
	sequences_out.open(options.sequences_out_file.c_str());
	substitutions_out.open(options.substitutions_out_file.c_str());

	substitutions_out << "Branch\tSubstitutions" << endl;
}

void Network::SampleParameters() {
	std::cout << "Sampling network parameters" << std::endl;
	SampleSequences();
}

/**
 * There are many different ancestral state sampling algorithms. This is the
 * simplest I could think of. Next, one might imagine allowing the ancestor to
 * be the same state as the descendant but choose the left or the right with
 * the probability of the substitution required. For example, then low
 * probability substitutions would not be required and the MCMC would sample
 * the high probability substitutions more often than the low probability
 * substitutions.
 * p(left) = p(substitution required if ancestor state is left state) / normalize
 * p(right) = p(substitution required if ancestor state is right state) / normalize
 */
void Network::SampleSequences() {
	if (not IsLeaf() and not is_constant) {
		for (int site = 0; site < sequence.size(); site++) {
			// < is correct because Random() returns a value between 0 and 0.99
			// inclusive.
			if (Random() < 0.5) {
				sequence.at(site) = left->sequence.at(site);
			} else {
				sequence.at(site) = right->sequence.at(site);
			}
		}
	}
}

/**
 * This is a simple sampling method to take the previous value into account.
 * The new value is allowed to be 50% greater or lower than the previous value.
 * This allows the new value to sample larger and smaller values.
 */
void Network::SampleDistances() {
	if (not is_constant)
		distance = distance * (0.5 + Random());
}

void Network::RecordState() {
	RecordSubNetworkState();
	AddGenerationEndIndicatorsToOutputFiles();
}

//STP: RecordState should be recursive so long as you traverse the Network post-
// order. Newick Network format is a post-order format. Order does not matter for
// the sequences and substitutions.
void Network::RecordSubNetworkState() {
	if (IsSubNetwork()) {
		Network_out << "(";
		left->RecordSubNetworkState();
		Network_out << ",";
		right->RecordSubNetworkState();
		Network_out << ")";
	}

	RecordNameAndDistanceToNetworkfile();
	RecordSequence();
	RecordSubstitutions();
}

/**
 * Need to remember to print Network at least once to attach the names of the
 * internal nodes to the Network.
 *
 * This needs a better name.
 */
//
void Network::RecordNameAndDistanceToNetworkfile() {
	Network_out << name << ":" << distance;
}

//STP: Recording the sequences requires decoding the sequences from integers to
// amino acids or nucleotides (chars).
void Network::RecordSequence() {
	sequences_out << ">" << name << endl;

	for (int site = 0; site < sequence.size(); site++) {
		sequences_out << integer_to_state[sequence.at(site)];
	}
	sequences_out << endl;
}

//STP: When are the substitutions sampled? Only when we sample the internal
// sequences. Then we could determine the substitutions once then. I don't
// want to figure out the substituitons again here.
void Network::RecordSubstitutions() {
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

	if (IsSubNetwork()) {
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
 * It might make reverseStateIndex a global since the Network needs to use it.
 * But only the Network needs to use it; no other part needs the information.
 * Should the Network be able to decode itself?
 *
 */
void Network::RecordChildSubstitutions(Network* child) {
	substitutions_out << child->name;
	for (int site = 0; site < sequence.size(); site++) {
		if (sequence.at(site) != child->sequence.at(site)) {
			//Note: This prints out one-indexed sequence positions.
			substitutions_out << "\t" << integer_to_state[sequence.at(site)]
					<< (site + 1) << integer_to_state[child->sequence.at(site)];
		}
	}
	substitutions_out << endl;
}

/**
 * Adding the final semicolon must happen at only the root node. But the
 * root node doesn't know it is the root node. So the model must call this
 * function after the root has recorded itself.
 *
 */
void Network::AddGenerationEndIndicatorsToOutputFiles() {
	Network_out << ";" << endl;
	sequences_out << "//" << endl;
	// This is required because the root node is not a child of anything and so
	// its substitutions are never printed. This would not be here if a Network
	// could record its own substitutions.
	substitutions_out << name << endl;
	substitutions_out << "//" << endl;
}


#include "Sequence.h"
#include <iostream>
#include "Environment.h"
#include <algorithm>
#include "IO.h"

extern double Random();
extern Environment env;
extern IO::Files files;

std::vector<std::string> aa({"A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"});
std::vector<std::string> nucleotides({"A", "T", "C", "G"});

std::ofstream SequenceAlignment::sequences_out;

// Sequence Alignment class.

SequenceAlignment::SequenceAlignment() {
	int states_option = env.get_int("states");

	switch(states_option) {
		case 0: // Nucleotide.
			states = nucleotides;
		case 1: // Amino acid.
			states = aa;
	}

	env.num_states = states.size();

	for(int i = 0; i < states.size(); i++) {
		state_to_integer[states[i]] = i;
		state_to_integer["-"] = -1;
		integer_to_state[i] = states[i];
		integer_to_state[-1] = "-";
	}

	env.state_to_integer = state_to_integer;
}

static const int gap_indicator = -1;

void SequenceAlignment::add(std::string name, std::string sequence_str) {
	// Adds extant sequence to alignment. This will not be sampled during the MCMC.
	std::vector<int> enc = EncodeSequence(sequence_str);
	taxa_names_to_sequences[name] = enc;
}

void SequenceAlignment::add(std::string name) {
	// Adds sequence to alignment, that WILL be sampled during MCMC. 
	// This is for the ancestral nodes.
	std::vector<int> enc = {};
	taxa_names_to_sequences[name] = enc;
}

void SequenceAlignment::print() {
	std::cout << "SEQUENCES" << std::endl;
	for(std::map<std::string, std::vector<int>>::iterator it = taxa_names_to_sequences.begin(); it != taxa_names_to_sequences.end(); ++it) {
		std::cout << ">" << it->first << "\n" << decodeSequence(it->second) << std::endl;
	}
}

void SequenceAlignment::Initialize() {
	// Process sequences.
	DetermineColumnsWithoutGaps();
	RemoveColumnsWithGapsFromSequences();

	// Setup output.
	files.add_file("sequences", env.get("sequences_out_file"), IOtype::OUTPUT);
	sequences_out = files.get_ofstream("sequences");

	// Setup Environment.
	env.n = (*taxa_names_to_sequences.begin()).second.size();
}

void SequenceAlignment::saveToFile(int gen, double l) {
	static int i = -1;
	++i;
	sequences_out << "#" << i << ":" << gen << ":" << l << std::endl;
	for(auto it = taxa_names_to_sequences.begin(); it != taxa_names_to_sequences.end(); ++it) {
		sequences_out << ">" << it->first << "\n" << decodeSequence(it->second) << std::endl;
	}
}

// Reading Fasta files.
std::vector<int> SequenceAlignment::EncodeSequence(const std::string &sequence) {
	/*
	 * Takes a string representation of a sequence and returns vector of integers.
	 * Also tracks the gaps in the alignment.
	 */
	std::vector<int> encoded_sequence(sequence.length());

	for (int site = 0; site < sequence.length(); site++) {
		std::string current_pos = sequence.substr(site, 1);
		encoded_sequence.at(site) = state_to_integer[current_pos];
	}
	return encoded_sequence;
}

void SequenceAlignment::DetermineColumnsWithoutGaps() {
	int number_of_sites = taxa_names_to_sequences.begin()->second.size();

	for (int site = 0; site < number_of_sites; site++) {
		// If site is not in columns with gaps
		if (columns_with_gaps.find(site) == columns_with_gaps.end()) {
			columns_without_gaps.push_back(site);
		}
	}
}

void SequenceAlignment::RemoveColumnsWithGapsFromSequences() {
	for (std::map<std::string, std::vector<int>>::iterator it = taxa_names_to_sequences.begin(); it != taxa_names_to_sequences.end(); ++it) {
		std::vector<int> encoded_sequence = it->second;
		it->second = RemoveGapsFromEncodedSequence(encoded_sequence);
	}
}

std::vector<int> SequenceAlignment::RemoveGapsFromEncodedSequence(std::vector<int> encoded_sequence) {
	std::vector<int> encoded_sequence_without_gaps(columns_without_gaps.size());

	for (int site = 0; site < columns_without_gaps.size(); site++) {
		encoded_sequence_without_gaps.at(site) = encoded_sequence.at(columns_without_gaps.at(site));
	}
	return encoded_sequence_without_gaps;
}

// Utilities
int SequenceAlignment::numCols() {
	return(columns_without_gaps.size());
}

std::string SequenceAlignment::decodeChar(int &c) {
	return(integer_to_state[c]);
}

std::string SequenceAlignment::decodeSequence(std::vector<int> &enc_seq) {
	std::string decoded_sequence;
	for(std::vector<int>::iterator it = enc_seq.begin(); it != enc_seq.end(); ++it) {
		decoded_sequence.append(decodeChar(*it));
	}
	return(decoded_sequence);
}

std::vector<int> SequenceAlignment::findParsimony(const std::vector<int> &s1, const std::vector<int> &s2) {
	std::vector<int> p = {};
	for(int i = 0; i < s1.size(); i++) {
		if(s1.at(i) == s2.at(i)) {
			//If states are the same.
			p.push_back(s1.at(i));
		} else {
			// If states are differant - check if one of the states is gap ("-").
			if(s1.at(i) == -1) {
				p.push_back(s2.at(i));
			} else if(s2.at(i) == -1) {
				p.push_back(s1.at(i));
			} else if(Random() < 0.5) {
				p.push_back(s1.at(i));
			} else {
				p.push_back(s2.at(i));
			}
		}
	}	
	return(p);
}

std::list<substitution> SequenceAlignment::findSubstitutions(const std::vector<int> &anc, const std::vector<int> &dec) {
	std::list<substitution> s = {};
	for(int i = 0; i < anc.size(); i++) {
		if(anc.at(i) != dec.at(i)) {
			substitution sub = {i, anc.at(i), dec.at(i)};
			s.push_back(sub);
		}
	}	
	return(s);
}


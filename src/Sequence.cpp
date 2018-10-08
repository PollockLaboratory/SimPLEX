#include "Sequence.h"
#include <iostream>
#include "Environment.h"

extern Environment env;

// Sequence Alignment class.

SequenceAlignment::SequenceAlignment() {
}

static const int gap_indicator = -1;

void SequenceAlignment::add(std::string name, std::string sequence_str) {
	std::cout << "Adding: " << name << std::endl;
	std::vector<int> enc = EncodeSequence(sequence_str);
	taxa_names_to_sequences[name] = new Sequence(enc, this);
}

std::vector<int> SequenceAlignment::EncodeSequence(std::string sequence) {
	/*
	 * Takes a string representation of a sequence and returns vector of integers.
	 * Also tracks the gaps in the alignment.
	 */
	std::cout << "Reading sequence: " << sequence << std::endl;
	std::vector<int> encoded_sequence(sequence.length());

	for (int site = 0; site < sequence.length(); site++) {
		std::string current_pos = sequence.substr(site, 1);
		if (sequence.find('-') != std::string::npos
				or sequence.find('?') != std::string::npos) { // When sequence.find == string::npos, the string was not found
			// There is gap.
			columns_with_gaps.insert(site);
			encoded_sequence.at(site) = gap_indicator;
		} else {
			//No gap.
			AddStateToStates(current_pos);
			encoded_sequence.at(site) = state_to_integer[current_pos];
		}
	}
	return encoded_sequence;
}

void SequenceAlignment::AddStateToStates(std::string state) {
	// Checks to see of state exists in state_to_interger_map.
	if (state_to_integer.find(state) == state_to_integer.end()) {
		int encoded_state = state_to_integer.size();
		state_to_integer[state] = encoded_state;
		integer_to_state[encoded_state] = state;

		states.push_back(state);
	}
}

void SequenceAlignment::print() {
	std::cout << "SEQUENCES" << std::endl;
	for(std::map<std::string, Sequence*>::iterator it = taxa_names_to_sequences.begin(); it != taxa_names_to_sequences.end(); ++it) {
		std::cout << "Name: " << it->first << " Sequence: " << (it->second)->as_str() << std::endl;
	}
}

void SequenceAlignment::DetermineColumnsWithoutGaps() {
	int number_of_sites = taxa_names_to_sequences.begin()->second->encoded_sequence.size();
	std::cout << "Number of columns in alignment: " << number_of_sites << std::endl;

	for (int site = 0; site < number_of_sites; site++) {
		// If site is not in columns with gaps
		if (columns_with_gaps.find(site) == columns_with_gaps.end()) {
			columns_without_gaps.push_back(site);
		}
	}

	std::cout << "Number of columns without gaps: " << columns_without_gaps.size() << std::endl;
}

void SequenceAlignment::RemoveColumnsWithGapsFromSequences() {
	if (env.debug) std::cout << "Removing gaps" << std::endl;

	for (std::map<std::string, Sequence*>::iterator it = taxa_names_to_sequences.begin(); it != taxa_names_to_sequences.end(); ++it) {
		std::vector<int> encoded_sequence = it->second->encoded_sequence;
		it->second->encoded_sequence = RemoveGapsFromEncodedSequence(encoded_sequence);
	}
}

std::vector<int> SequenceAlignment::RemoveGapsFromEncodedSequence(std::vector<int> encoded_sequence) {
	std::vector<int> encoded_sequence_without_gaps(columns_without_gaps.size());

	for (int site = 0; site < columns_without_gaps.size(); site++) {
		encoded_sequence_without_gaps.at(site) = encoded_sequence.at(columns_without_gaps.at(site));
	}
	return encoded_sequence_without_gaps;
}

void SequenceAlignment::Initialize() {
	DetermineColumnsWithoutGaps();
	RemoveColumnsWithGapsFromSequences();
}

// Indervidual Sequence class.

Sequence::Sequence(std::vector<int> enc, SequenceAlignment* MSA) {
	encoded_sequence = enc;
	this->MSA = MSA;
}

std::string Sequence::as_str() {
	std::string decoded_sequence;
	for(std::vector<int>::iterator it = encoded_sequence.begin(); it != encoded_sequence.end(); ++it) {
		decoded_sequence.append(MSA->integer_to_state[*it]);
	}
	return(decoded_sequence);
}



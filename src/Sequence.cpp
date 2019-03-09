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

SequenceAlignment::SequenceAlignment(const States* states) {
  env.num_states = states->n;
  state_to_integer = states->state_to_int;
  integer_to_state = states->int_to_state;
  env.state_to_integer = state_to_integer;
}

SequenceAlignment::SequenceAlignment(const SequenceAlignment &msa) {
  // The copy constructor.
  // Check whether this is truely copying.
  taxa_names_to_sequences = msa.taxa_names_to_sequences;
  columns_with_gaps = msa.columns_with_gaps;
  columns_without_gaps = msa.columns_without_gaps;
  states = msa.states;
  state_to_integer = msa.state_to_integer;
  integer_to_state = msa.integer_to_state;
  MSA_list = msa.MSA_list;
  current_MSA = MSA_list->begin();
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

void SequenceAlignment::Initialize(std::list<SequenceAlignment*>* msa_List) {
  // Process sequences.
  DetermineColumnsWithoutGaps();
  RemoveColumnsWithGapsFromSequences();
  MSA_list = msa_List;
  // current_MSA = MSA_list->begin();

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

  for (unsigned int site = 0; site < sequence.length(); site++) {
    std::string current_pos = sequence.substr(site, 1);
    encoded_sequence.at(site) = state_to_integer[current_pos];
  }
  return(encoded_sequence);
}

void SequenceAlignment::DetermineColumnsWithoutGaps() {
  unsigned int number_of_sites = taxa_names_to_sequences.begin()->second.size();

  for(unsigned int site = 0; site < number_of_sites; site++) {
    // If site is not in columns with gaps
    if(columns_with_gaps.find(site) == columns_with_gaps.end()) {
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

  for (unsigned int site = 0; site < columns_without_gaps.size(); site++) {
    encoded_sequence_without_gaps.at(site) = encoded_sequence.at(columns_without_gaps.at(site));
  }
  return(encoded_sequence_without_gaps);
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
  for(unsigned int i = 0; i < s1.size(); i++) {
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

//std::list<substitution> SequenceAlignment::findSubstitutions(const std::vector<int> &anc, const std::vector<int> &dec) {
// std::list<substitution> s = {};
//  for(int i = 0; i < anc.size(); i++) {
//   if(anc.at(i) != dec.at(i)) {
//     substitution sub = {i, anc.at(i), dec.at(i)};
//     s.push_back(sub);
//   }
// }
// return(s);
//}

std::list<std::string> SequenceAlignment::getNodeNames() {
  std::list<std::string> names;
  for(auto it = taxa_names_to_sequences.begin(); it != taxa_names_to_sequences.end(); ++it) {
    names.push_back(it->first);
  }
  return(names);
}

// Ancestral Sequences
void SequenceAlignment::step_to_next_MSA() {
  if(env.ancestral_sequences != 1) {
    std::cerr << "Error: Invalid step_to_next_MSA call, ancestral sequences not previously determined." << std::endl;
    exit(EXIT_FAILURE);
  }

  ++current_MSA;
  if(current_MSA == MSA_list->end()) {
    current_MSA = MSA_list->begin();
  }

  for(auto it = taxa_names_to_sequences.begin(); it != taxa_names_to_sequences.end(); ++it) {
    // std::cout << "Name: " << it->first << std::endl;
    // std::cout << "Current sequence: " << decodeSequence(it->second) << std::endl;
    // std::cout << "New Sequence: " << decodeSequence((*current_MSA)->taxa_names_to_sequences[it->first]) << std::endl;
    it->second = (*current_MSA)->taxa_names_to_sequences[it->first];
  }
}


#ifndef Sequence_h_
#define Sequence_h_

#include<string>
#include<map>
#include<vector>
#include<list>
#include<set>

#include "SubstitutionModels/SubstitutionModel.h"
#include "../IO/SequencesParser.h"

struct substitution {
	int pos;
	int anc;
	int dec;
};

class SequenceAlignment {
 public:
  SequenceAlignment(const States*);
  SequenceAlignment(const SequenceAlignment &msa);
  std::map<std::string, std::vector<int>> taxa_names_to_sequences;

  std::set<int> columns_with_gaps;
  std::vector<int> columns_without_gaps;

  std::vector<std::string> states;
  std::map<std::string, int> state_to_integer;
  std::map<int, std::string> integer_to_state;

  // Known ancestral sequences.
  std::list<SequenceAlignment*> *MSA_list; // List of other sequence alignments that could be used.
  std::list<SequenceAlignment*>::iterator current_MSA; // Iterator pointer to the MSA in use.
  void step_to_next_MSA();

  // Adding sequences to alignment.
  void add(std::string name, std::string sequence_str);
  void add(std::string name);

  void print();
  void Initialize(std::list<SequenceAlignment*>*);
  void Initialize(IO::RawMSA* &raw_msa);
  void saveToFile(int gen, double l);

  // Utilities
  int numCols();
  std::string decodeChar(int &c);
  std::string decodeSequence(std::vector<int> &enc_seq);
  static std::vector<int> findParsimony(const std::vector<int> &s1, const std::vector<int> &s2);
  std::list<std::string> getNodeNames();
 private:
  // Processing input sequences - fasta.
  std::vector<int> EncodeSequence(const std::string &sequence);
  void DetermineColumnsWithoutGaps();
  void RemoveColumnsWithGapsFromSequences();
  std::vector<int> RemoveGapsFromEncodedSequence(std::vector<int> encoded_sequence);

  // Output.
  static std::ofstream sequences_out;
};
#endif

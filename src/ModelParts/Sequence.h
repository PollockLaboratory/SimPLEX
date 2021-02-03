#ifndef Sequence_h_
#define Sequence_h_

#include <string>
#include <map>
#include <vector>
#include <list>
#include <set>

#include "SubstitutionModels/SubstitutionModel.h"
#include "../IO/SequencesParser.h"

// Temporary
//class Tree;

class TreeNode;

struct substitution {
	int pos;
	int anc;
	int dec;
};

class SequenceAlignment {
 public:
  Tree* tree;
  unsigned int n_columns;
  unsigned int n_states;
  SequenceAlignment(const States*);
  SequenceAlignment(const SequenceAlignment &msa);

  std::map<std::string, std::vector<int>> taxa_names_to_sequences;
  std::map<std::string, std::vector<bool>> taxa_names_to_gaps;
  std::list<std::string> base_sequences;

  std::set<std::string> states;
  std::map<std::string, int> state_to_integer;
  std::map<int, std::string> integer_to_state;

  // Known ancestral sequences.
  std::list<SequenceAlignment*> *MSA_list; // List of other sequence alignments that could be used.
  std::list<SequenceAlignment*>::iterator current_MSA; // Iterator pointer to the MSA in use.
  void step_to_next_MSA();

  // Adding sequences to alignment.
  void add(std::string name, std::string sequence_str);
  void add(std::string name);
  void add_base(std::string name, std::string sequence_str);

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

  // New
  void syncWithTree(Tree* tree);
  void identify_gaps();
  sample_status sample();
 private:
  // Processing input sequences - fasta.
  std::vector<int> EncodeSequence(const std::string &sequence);

  // Sampling - new
  std::map<std::string, float**> taxa_names_to_state_probs;

  void sample_base_sequences(std::list<int>);
  void calculate_state_probabilities(TreeNode*, std::list<int>);
  void calculate_state_probabilities_pos(TreeNode*, unsigned int, TreeNode*, TreeNode*, TreeNode*);
  int pick_state_from_probabilities(TreeNode*, int);
};

class SequenceAlignmentParameter : public SampleableComponent {
private:
  SequenceAlignment* msa;
public:
  SequenceAlignmentParameter(SequenceAlignment*);

  void print() override;
  std::string get_type() override;
  
  sample_status sample() override;
  void undo() override;
  void fix() override;
  void refresh() override;

  std::string get_state_header() override;
  std::string get_state() override;

  void save_to_file(int gen, double l);
};

#endif

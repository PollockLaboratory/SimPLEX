#ifndef Sequence_h_
#define Sequence_h_

#include <string>
#include <map>
#include <vector>
#include <list>
#include <set>
#include <boost/multiprecision/cpp_int.hpp>

#include "AbstractComponent.h"
#include "../IO/SequencesParser.h"
#include "SubstitutionModels/States.h"

using boost::multiprecision::uint128_t;

class Tree;
class TreeNode;
class BranchSegment;

class SequenceAlignment {
 public:
  Tree* tree;
  std::string domain_name; // Domain.
  unsigned int n_columns;
  unsigned int n_states;
  SequenceAlignment(std::string name, std::string msa_out, std::string subs_out, const States*);

  // States
  std::set<std::string> states;
  std::map<std::string, state_element> state_element_encode;
  std::map<state_element, std::string> state_element_decode;

  // Adding sequences to alignment.
  void add(std::string name, std::string sequence_str);
  void add(std::string name);
  void add_base(std::string name, const IO::FreqSequence &seq);

  void print();
  void Initialize(IO::RawMSA raw_msa);
  void saveToFile(int i, uint128_t gen, double l);

  // Utilities
  int n_cols();
  std::string decode_state_element(state_element c);
  std::string decode_state_element_sequence(std::vector<state_element> &enc_seq);
  static std::vector<signed char> find_parsimony(const std::vector<state_element> &s1, const std::vector<state_element> &s2);
  std::list<std::string> getNodeNames();

  // Validation.
  bool match_structure(SequenceAlignment*);
  bool validate(std::list<std::string> seq_names, std::map<std::string, SequenceAlignment*> other_alignments);

  void syncWithTree(std::string name, unsigned int id, Tree* tree);
  void identify_gaps();

  sample_status sample(const std::list<unsigned int>&);
 private:
  // Processing input sequences - fasta.
  std::vector<signed char> encode_sequence(const std::string &sequence);

  // Probability distributions of each state.
  // Taxa name -> residue position -> state -> probability.
  // double** is a matrix where i is the sequence position and j is the state id.
  std::map<std::string, double**> prior_state_distribution; // These are the priors.
  std::map<std::string, double**> marginal_state_distribution;

  // Sequences
  std::map<std::string, std::vector<state_element>> taxa_names_to_sequences;

  // Gaps
  std::map<std::string, std::vector<bool>> taxa_names_to_gaps;

  void reset_to_base(std::string name, const std::list<unsigned int>& positions);
  void normalize_state_probs(TreeNode* node, unsigned int pos);

  // Marginal Calculations for indervidual positions.
  double find_state_prob_given_dec_branch(BranchSegment* branch, state_element state_i, double* state_probs, std::vector<Valuable*> rv, double u, unsigned int pos);
  double find_state_prob_given_anc_branch(BranchSegment* branch, state_element state_i, double* state_probs, TreeNode* node, double u, unsigned int pos);
  void find_marginal_at_pos(TreeNode*, unsigned int, TreeNode*, TreeNode*, TreeNode*);

  // Marginal posterior calculations for whole sequences.
  void find_state_probs_dec_only(TreeNode*, std::list<unsigned int>); // First Recursion.
  void update_state_probs(TreeNode* node, unsigned int pos, TreeNode* up_node); // Second Recursion
  void find_state_probs_all(TreeNode*, std::list<unsigned int>); // Third Recursion

  // Optimize
  void fast_update_state_probs_tips(TreeNode* node, unsigned int pos, TreeNode* up_node); // Second Recursion

  // Picking states
  int pick_state_from_probabilities(TreeNode*, int);
  void pick_states_for_node(TreeNode*, const std::list<unsigned int>&);

  void reconstruct_expand(const std::list<TreeNode*>&, const std::list<unsigned int>&);

  // Outputs
  std::string seqs_out_identifier;
  std::string seqs_out_file;
  std::string substitutions_out_identifier;
  std::string substitutions_out_file;
};

class SequenceAlignmentParameter : public SampleableComponent {
private:
  SequenceAlignment* msa;
  unsigned int n_sample; // Number of positions sampled each time.
  unsigned int n_cols; // Number of columns in the MSA.
  unsigned int sample_loc; // position of last sample.

  int save_count;
public:
  SequenceAlignmentParameter(SequenceAlignment*, unsigned int n_sample);

  void print() override;
  std::string get_type() override;
  
  sample_status sample() override;
  void undo() override;
  void fix() override;
  void refresh() override;

  std::string get_state_header() override;
  std::string get_state() override;

  void save_to_file(uint128_t gen, double l);
};

#endif

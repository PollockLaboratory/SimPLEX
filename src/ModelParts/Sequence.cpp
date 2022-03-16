#include <sstream>
#include <algorithm>

#include "Sequence.h"
#include "../Environment.h"
#include "../IO/Files.h"

#include "Trees/TreeParts.h"
#include "Trees/Tree.h"

// Globals
extern double Random();
extern Environment env;
extern IO::Files files;

std::vector<std::string> aa({"A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"});
std::vector<std::string> nucleotides({"A", "T", "C", "G"});

static const signed char gap_indicator = -1;

// Sequence Alignment class.

SequenceAlignment::SequenceAlignment(std::string name, std::string msa_out, std::string subs_out, const States* states) : domain_name(name) {
  this->states = states->possible;
  n_states = states->n;

  state_to_integer = states->state_to_int;
  integer_to_state = states->int_to_state;

  seqs_out_file = msa_out;
  substitutions_out_file = subs_out;
}

void SequenceAlignment::add(std::string name, std::string sequence_str) {
  // Adds extant sequence to alignment. This will not be sampled during the MCMC.
  std::vector<signed char> enc = encode_sequence(sequence_str);
  taxa_names_to_sequences[name] = enc;
}

void SequenceAlignment::add(std::string name) {
  // Adds sequence to alignment, that WILL be sampled during MCMC.
  // This is for the ancestral nodes.
  std::vector<signed char> enc = {};
  taxa_names_to_sequences[name] = enc;
}

float** create_state_probability_vector(unsigned int n_cols, unsigned int n_states) {
  float** m = new float*[n_cols];
  for(unsigned int i = 0; i < n_cols; i++) {
    m[i] = new float[n_states];
    for(unsigned int j = 0; j < n_states; j++) {
      m[i][j] = 0.0;
    }
  }

  return(m);
}

void SequenceAlignment::add_base(std::string name, const IO::FreqSequence &seq) {
  base_sequences.push_front(name);
  add(name, sequenceAsStr_highestFreq(seq));

  base_taxa_state_probs[name] = create_state_probability_vector(seq.size(), n_states);

  int pos = 0;
  for(auto it = seq.begin(); it != seq.end(); ++it) {
    for(auto jt = it->begin(); jt != it->end(); ++jt) {
      if(jt->state != '-') {
	base_taxa_state_probs[name][pos][state_to_integer[std::string(1, jt->state)]] = jt->freq;
      }
    }
    pos++;
  }
}

void SequenceAlignment::print() {
  std::cout << "SEQUENCES" << std::endl;
  for(std::map<std::string, std::vector<signed char>>::iterator it = taxa_names_to_sequences.begin(); it != taxa_names_to_sequences.end(); ++it) {
    std::cout << ">" << it->first << "\n" << decodeSequence(it->second) << std::endl;
  }
}

void SequenceAlignment::Initialize(IO::RawMSA raw_msa) {
  for(auto it = raw_msa.seqs.begin(); it != raw_msa.seqs.end(); ++it) {
    add_base(it->first, it->second);
  }

  // Setup output.
  seqs_out_identifier = domain_name + "_sequences_out";
  files.add_file(seqs_out_identifier, seqs_out_file, IOtype::OUTPUT);

  substitutions_out_identifier = domain_name + "_substitutions_out";
  files.add_file(substitutions_out_identifier, substitutions_out_file, IOtype::OUTPUT);
  files.write_to_file(substitutions_out_identifier, "I,GEN,LogL,Ancestral,Decendant,Substitutions\n");

  n_columns = (*taxa_names_to_sequences.begin()).second.size();
}

void SequenceAlignment::saveToFile(int save_count, uint128_t gen, double l) {
  std::ostringstream buffer;
  buffer << "#" << save_count << ":" << gen << ":" << l << std::endl;
  for(auto it = taxa_names_to_sequences.begin(); it != taxa_names_to_sequences.end(); ++it) {
    buffer << ">" << it->first << "\n" << decodeSequence(it->second) << std::endl;
  }

  files.write_to_file(seqs_out_identifier, buffer.str());

  std::ostringstream subs_buffer;

  // Save the substitutions for a particular state.
  std::list<BranchSegment*> branches = tree->get_branches();
  for(auto it = branches.begin(); it != branches.end(); ++it) {
    subs_buffer << save_count << "," << gen << "," << l << ",";
    subs_buffer << (*it)->ancestral->name << "," << (*it)->decendant->name << ",[ ";
    std::vector<Substitution> subs = (*it)->get_substitutions(domain_name);
    for(unsigned int pos = 0; pos < subs.size(); pos++) {
      if(subs[pos].occuredp == true) {
	int anc = (*it)->ancestral->sequences[domain_name]->at(pos);
	int dec = (*it)->decendant->sequences[domain_name]->at(pos);

	// Includes virtual substitutions.
	subs_buffer << integer_to_state[anc] << pos << integer_to_state[dec] << " ";
      }
    }
    subs_buffer << "]\n";
  }
  files.write_to_file(substitutions_out_identifier, subs_buffer.str());
}

void SequenceAlignment::syncWithTree(std::string name, unsigned int id, Tree* tree) {
  /*
    Connects all the tree nodes on the matching sequences in the MSAs for each state domain.
    Add new sequences to MSA for nodes present on the tree but missing in the alignment.
    This results in pointers from the TreeNodes to the sequence vectors.
   */

  this->tree = tree;
  std::cout << "\t\tSyncing " << name << " states to tree. ID: " << id << std::endl;
  std::list<TreeNode*> nodes = tree->nodes();

  for(auto it = nodes.begin(); it != nodes.end(); ++it) {
    TreeNode* n = *it;

    if(taxa_names_to_sequences.count(n->name)) {
      n->sequences[name] = &(taxa_names_to_sequences.at(n->name));
    } else {
      if(n->isTip()){
	std::cerr << "Error: Missing sequence for \"" << n->name << "\"." << std::endl;
	exit(EXIT_FAILURE);
      } else {
	// Add new sequence to sequence alignments.
	add(n->name);
	n->sequences[name] = &(taxa_names_to_sequences.at(n->name));
	  
	// Fill missing sequences/
	if(n->left != 0 and n->right == 0) {
	  // Internal Continous.
	  TreeNode* dsNode = n->left->decendant; // ds = downstream.
	  *(n->sequences[name]) = *(dsNode->sequences[name]);
	} else {
	  // Root or internal branch.
	  vector<signed char> dsNodeLseq = *(n->left->decendant->sequences[name]);
	  vector<signed char> dsNodeRseq = *(n->right->decendant->sequences[name]);
	  vector<signed char> p = findParsimony(dsNodeLseq, dsNodeRseq);
	  *(n->sequences[name]) = p;
	}
      }
    }
  }

  identify_gaps();
}
					     
void SequenceAlignment::identify_gaps() {
  std::list<TreeNode*> nodes = tree->nodes();
  for(auto node = nodes.begin(); node != nodes.end(); ++node) {
    //std::cout << (*node)->name << std::endl;

    std::vector<signed char> sequence = taxa_names_to_sequences[(*node)->name];
    std::vector<bool> gaps(sequence.size(), false);

    taxa_names_to_state_probs[(*node)->name] = create_state_probability_vector(n_columns, n_states);

    if((*node)->isTip()) {
      for(unsigned int pos = 0; pos < sequence.size(); pos++) {
	if(sequence.at(pos) == -1) {
	  gaps[pos] = true;
	} else {
	  gaps[pos] = false;
	}
      }
    } else {
      for(unsigned int pos = 0; pos < sequence.size(); pos++) {
	// Checks left branch first - there is always a left branch, except tips.
	if(taxa_names_to_gaps[(*node)->left->decendant->name][pos]) {
	  // Checks right branch next - right branches only on branching segment.
	  if((*node)->right) {
	    if(taxa_names_to_gaps[(*node)->right->decendant->name][pos]) {
	      gaps[pos] = true;
	    } else {
	      gaps[pos] = false;
	    }
	  } else {
	    gaps[pos] = true;
	  }
	} else {
	  gaps[pos] = false;
	}
      }
    }
    taxa_names_to_gaps[(*node)->name] = gaps;
  }
}

// Reading Fasta files.
std::vector<signed char> SequenceAlignment::encode_sequence(const std::string &sequence) {
  /*
   * Takes a string representation of a sequence and returns vector of integers.
   * Also tracks the gaps in the alignment.
   */
  std::vector<signed char> encoded_sequence(sequence.length());

  for (unsigned int site = 0; site < sequence.length(); site++) {
    std::string current_pos = sequence.substr(site, 1);
    try {
      encoded_sequence.at(site) = state_to_integer.at(current_pos);
    } catch(const std::out_of_range& e) {
      // Temporary solution.
      if(current_pos == "-"){
	encoded_sequence.at(site) = -1;
      } else {
	std::cerr << "Error: state \"" << current_pos << "\" in sequence alignment is not recognised. " << std::endl;
	exit(EXIT_FAILURE);
      }
    }
  }
  
  return(encoded_sequence);
}

// Utilities
int SequenceAlignment::numCols() {
  return(n_columns);
}

std::string SequenceAlignment::decodeChar(signed char c) {
  return(integer_to_state[c]);
}

std::string SequenceAlignment::decodeSequence(std::vector<signed char> &enc_seq) {
  std::string decoded_sequence;
  for(std::vector<signed char>::iterator it = enc_seq.begin(); it != enc_seq.end(); ++it) {
    decoded_sequence.append(decodeChar(*it));
  }
  return(decoded_sequence);
}

std::vector<signed char> SequenceAlignment::findParsimony(const std::vector<signed char> &s1, const std::vector<signed char> &s2) {
  std::vector<signed char> p = {};
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

std::list<std::string> SequenceAlignment::getNodeNames() {
  std::list<std::string> names;
  for(auto it = taxa_names_to_sequences.begin(); it != taxa_names_to_sequences.end(); ++it) {
    names.push_back(it->first);
  }
  return(names);
}

void SequenceAlignment::reset_base_probabilities() {
  for(auto seq_name = base_sequences.begin(); seq_name != base_sequences.end(); ++seq_name) {
    for(unsigned int i = 0; i < n_columns; i++) {
      for(unsigned int j = 0; j < n_states; j++) {
	taxa_names_to_state_probs[*seq_name][i][j] = base_taxa_state_probs[*seq_name][i][j];
      }
    }
  }
}

void SequenceAlignment::normalize_state_probs(TreeNode* node, unsigned int pos) {
  double normalize_total = 0.0;

  for(unsigned int i = 0; i < n_states; i++) {
    normalize_total += taxa_names_to_state_probs[node->name][pos][i];
  }

  if(normalize_total != 0.0) {
    for(unsigned int i = 0; i < n_states; i++) {
      taxa_names_to_state_probs[node->name][pos][i] /= normalize_total;
    }
  }
}

// Calculating Probabilities.
inline double calc_substitution_prob(double rate, float t_b, double u) {
  /*
   * rate = substitution rate.
   * t_b = branch length.
   * u = uniformisation constant.
   */
  return((rate * t_b) / (1.0 + (u * t_b)));
}

inline double calc_no_substitution_prob(double rate, float t_b, double u) {
  // Probability that there is any virtual substitution | given no substitution - equation (9) Rapid Likelihood Analysis on Large Phylogenies.
  double prob_virtual = 1.0 - (1.0 / (1.0 + (rate * t_b)));
  double denom = 1.0 / (1.0 + (u * t_b));

  //     Virtual Substitution                       No substitution
  return((prob_virtual * ((rate * t_b) * denom)) + ((1.0 - prob_virtual) * denom));
}

float SequenceAlignment::find_state_prob_given_dec_branch(unsigned char state_i, float* state_probs, std::vector<Valuable*> rv, float t_b, double u) {
  /*
   * Find state probability given decendent branch
   * state_probs = the marginal posterior distribution of the state at the node below.
   */
  float prob = 0.0;

  for(signed char state_j = 0; state_j < (signed char)n_states; state_j++) {
    double rate = rv[state_j]->get_value();
    double state_prob = state_probs[state_j];
    
    if(state_i != state_j) {
      // Normal Substitution
      prob += state_prob * calc_substitution_prob(rate, t_b, u); // Probability of the substitution.
    } else {
      // No substitution - possibly virtual.
      prob += state_prob * calc_no_substitution_prob(rate, t_b, u);
    }
  }

  return(prob);
}

float SequenceAlignment::find_state_prob_given_anc_branch(unsigned char state_j, float* state_probs, TreeNode* up_node, float t_b, double u, unsigned int pos) {
  /*
   * Find state probability given ancestor branch.
   */

  float prob = 0.0;

  for(signed char state_i = 0; state_i < (signed char)n_states; state_i++) {
    unsigned long extended_state = up_node->get_hypothetical_hash_state(pos, domain_name, state_i);
    RateVector* rv = up_node->SM->selectRateVector({pos, domain_name, extended_state});

    // Is this the right rate.
    double rate = rv->rates[state_j]->get_value(); // i -> j rate.
    double state_prob = state_probs[state_i];

    if(state_i != state_j) {
      // Normal Substitution.
      prob += state_prob * calc_substitution_prob(rate, t_b, u); // Probability of the substitution.
    } else {
      // No substition - possibly virtual.
      prob += state_prob * calc_no_substitution_prob(rate, t_b, u);
    }
  }

  return(prob);
}

void SequenceAlignment::find_new_state_probs_at_pos(TreeNode* node, unsigned int pos, TreeNode* left_node, TreeNode* right_node, TreeNode* up_node) {
  double u = node->SM->get_u();
  RateVector* rv;
  float left_prob = 1.0;
  float right_prob = 1.0;
  float up_prob = 1.0;

  for(signed char state_i = 0; state_i < (signed char)n_states; state_i++) {
    unsigned long extended_state = node->get_hypothetical_hash_state(pos, domain_name, state_i);

    rv = node->SM->selectRateVector({pos, domain_name, extended_state});

    // Contribution of left branch.
    if(left_node != nullptr) {
      if(not taxa_names_to_gaps[left_node->name].at(pos)) {
	left_prob = find_state_prob_given_dec_branch(state_i, taxa_names_to_state_probs[left_node->name][pos], rv->rates, left_node->distance, u);
      } else {
	  left_prob = 1.0;
      }
    } else {
      left_prob = 1.0;
    }

    // Contribution of right branch.
    if(right_node != nullptr) {
      if(not taxa_names_to_gaps[right_node->name][pos]) {
	right_prob = find_state_prob_given_dec_branch(state_i, taxa_names_to_state_probs[right_node->name][pos], rv->rates, right_node->distance, u);
      } else {
	right_prob = 1.0;
      } 
    } else {
      right_prob = 1.0;
    }

    if(up_node != nullptr) {
      if(not taxa_names_to_gaps[up_node->name][pos]) {
	up_prob = find_state_prob_given_anc_branch(state_i, taxa_names_to_state_probs[up_node->name][pos], up_node, up_node->distance, u, pos);
      }
    } else {
      up_prob = 1.0;
    }

    float total = left_prob * right_prob * up_prob;
    taxa_names_to_state_probs[node->name][pos][state_i] = total;
  }
}

void SequenceAlignment::calculate_state_probabilities(TreeNode* node, std::list<unsigned int> positions) {
  /*
   * Finds the marginal posterior distribution for each position at a given node.
   * Only uses infomation from nodes below - used for upward recursion.
   */
  
  std::string name = node->name;
  std::vector<bool> gaps = taxa_names_to_gaps[name];
  if(not node->isTip()) {
    // Node may or may not be a branch node - therefore may only have one child which is always the left one.
    // Set to nullptr if no right branch;
    TreeNode* right_node;
    if(node->right) {
      right_node = node->right->decendant;
    } else {
      right_node = nullptr;
    }

    for(auto pos = positions.begin(); pos != positions.end(); ++pos) {
      if(not gaps[*pos]) {
	find_new_state_probs_at_pos(node, *pos, node->left->decendant, right_node, nullptr);
	normalize_state_probs(node, *pos);
      }
    }
  }
}

// TODO refactor.
void SequenceAlignment::calculate_state_probabilities_full(TreeNode* node, std::list<unsigned int> positions) {
  /*
   * Finds the marginal posterior distribution for each position at a given node.
   * Only uses infomation from nodes below - used for upward recursion.
   */
  
  std::string name = node->name;
  std::vector<bool> gaps = taxa_names_to_gaps[name];
  if(not node->isTip()) {
    // Node may or may not be a branch node - therefore may only have one child which is always the left one.
    // Set to nullptr if no right branch;
    TreeNode* right_node;
    if(node->right) {
      right_node = node->right->decendant;
    } else {
      right_node = nullptr;
    }

    TreeNode* up_node;
    if(node->up) {
      up_node = node->up->ancestral;
    } else {
      up_node = nullptr;
    }

    for(auto pos = positions.begin(); pos != positions.end(); ++pos) {
      if(not gaps[*pos]) {
	find_new_state_probs_at_pos(node, *pos, node->left->decendant, right_node, up_node);
	normalize_state_probs(node, *pos);
      }
    }
  }
}

void SequenceAlignment::update_state_probs(TreeNode* node, unsigned int pos, TreeNode* up_node) {
  /*
   *
   */

  double u = node->SM->get_u();
  float up_prob = 0.0;

  for(signed char state_j = 0; state_j < (signed char)n_states; state_j++) {
    if(up_node != nullptr) {
      float t_b = node->distance;

      up_prob = find_state_prob_given_anc_branch(state_j, taxa_names_to_state_probs[up_node->name][pos], up_node, t_b, u, pos);
    } else {
      up_prob = 1.0;
    }

    taxa_names_to_state_probs[node->name][pos][state_j] *= up_prob;
  }
}

int SequenceAlignment::pick_state_from_probabilities(TreeNode* node, int pos) {
  /*
   * Picks a state from the marginal posterior distribution (taxa_names _to_state_probs).
   */
  float* probs = taxa_names_to_state_probs[node->name][pos];
 
  double r = Random();
  double acc = 0.0;
  int val = -1;
  for(unsigned int i = 0; i < n_states; i++) {
    acc += probs[i];
    if(r < acc and val == -1) {
      val = i;
      probs[i] = 1.0;
    } else {
      probs[i] = 0.0;
    }
  }

  if(val == -1) {
    std::cerr << "Error: incorrectly picking a state: " <<  pos << " " << node->name << std::endl;
    exit(EXIT_FAILURE);
  }

  return(val);
}

sample_status SequenceAlignment::sample() {
  const std::list<TreeNode*> nodes = tree->nodes();

  // Makes list of all positions.
  // Leaving room for a feature where not all positions are sampled each time.
  std::list<unsigned int> positions = {};
  for(unsigned int i = 0; i < n_columns; i++) {
    positions.push_back(i);
  };

  // This is important as states at tips can be uncertain.
  reset_base_probabilities();

  // Find state probabilities.
  // Nodes are ordered in the list such that they are visted in order up the tree.
  // Upward recursion.
  for(auto n = nodes.begin(); n != nodes.end(); ++n) {
    if(not (*n)->isTip()) {
      calculate_state_probabilities(*n, positions);
    }
  }

  // Reverse recursion.
  // Skip first element of reverse list as thats the root - no need to sample second time.
  for(auto n = nodes.rbegin(); n != nodes.rend(); ++n) {
    TreeNode* node = *n;
    std::vector<bool> gaps = taxa_names_to_gaps[node->name];

    // Reculaculate state probability vector - including up branch.
    TreeNode* up_node;
    if(node->up) {
      up_node = node->up->ancestral;
    } else {
      up_node = nullptr;
    }

    for(auto pos = positions.begin(); pos != positions.end(); ++pos) {
      if((not gaps[*pos]) and (not (up_node == nullptr))) {
	update_state_probs(node, *pos, up_node);
	normalize_state_probs(node, *pos);
      }
    }
  }

  // 3rd Recursion.
  for(auto n = nodes.rbegin(); n != nodes.rend(); ++n) {
    TreeNode* node = *n;
    std::vector<bool> gaps = taxa_names_to_gaps[node->name];

    for(auto pos = positions.begin(); pos != positions.end(); ++pos) {
      if(not (*n)->isTip()) {
	calculate_state_probabilities(*n, positions);
      }

      // Pick sequence.
      if(gaps[*pos]) {
	taxa_names_to_sequences[(*n)->name][*pos] = -1;
      } else {
	taxa_names_to_sequences[(*n)->name][*pos] = pick_state_from_probabilities(node, *pos);
      }
    } 
  }

  // Test
  //TreeNode* node = nodes[rand() % nodes.size()];
  //std::cout << node->name << std::endl;

  return(sample_status({false, true, true}));
}

// This need to be much more robust - add better error handling.
bool SequenceAlignment::validate(std::list<std::string> seq_names, std::map<std::string, SequenceAlignment*> other_alignments) {
  for(auto n = seq_names.begin(); n != seq_names.end(); ++n) {
    if(taxa_names_to_sequences.find(*n) == taxa_names_to_sequences.end()) {
      std::cerr << "Error: sequence alignment " << domain_name << " is missing sequence for " << *n << std::endl;
      exit(EXIT_FAILURE);
    }
  }
  return(true);
}

// This is not used currently.
bool SequenceAlignment::match_structure(SequenceAlignment* cmp_msa) {
  for(auto seq = taxa_names_to_sequences.begin(); seq != taxa_names_to_sequences.end(); ++seq) {
    // Check there are corresponding sequences.
    auto it = cmp_msa->taxa_names_to_sequences.find(seq->first);
    if(it != cmp_msa->taxa_names_to_sequences.end()) {
      //Check length of sequence matches.
      if(seq->second.size() != it->second.size()) {
	std::cerr << "Error: in sequence alignment \"" << domain_name << "\": sequences are not the same length as reference." << std::endl;
	exit(EXIT_FAILURE);
	return(false);
      } else {
	for(unsigned int i = 0; i < seq->second.size(); i++) {
	  if((seq->second[i] == -1) xor (it->second[i] == -1)) {
	    std::cout << seq->second[i] << " " << it->second[i] << std::endl;
	    std::cerr << "Error: in sequence alignment \"" << domain_name << "\" in sequence " << seq->first << " at position " << i << " inconsistant gaps." << std::endl;
	    exit(EXIT_FAILURE);
	    return(false);
	  }
	} 
      }
    } else {
      std::cerr << "Error: in sequence alignment \"" << domain_name << "\": sequence for \"" << seq->first << "\" is not found in reference." << std::endl;
      return(false);
    }
  }
  return(true);
}

SequenceAlignmentParameter::SequenceAlignmentParameter(SequenceAlignment* msa) : SampleableComponent("SequenceAlignment") {
  save_count = -1;
  this->msa = msa;
}

void SequenceAlignmentParameter::print() {
  std::cout << "SequenceAlignment-" << msa->domain_name << std::endl;
}

std::string SequenceAlignmentParameter::get_type() {
  return("SEQUENCE_ALIGNMENT");
}

sample_status SequenceAlignmentParameter::sample() {
  std::cout << "Sampling: " << msa->domain_name << std::endl;
  sample_status s = msa->sample();
  return(s);
}

void SequenceAlignmentParameter::undo() {
  std::cerr << "Error: SequenceAlignmentSampling cannot be undone." << std::endl;
  exit(EXIT_FAILURE);
}

void SequenceAlignmentParameter::fix() {
}

void SequenceAlignmentParameter::refresh() {
}

std::string SequenceAlignmentParameter::get_state_header() {
  return(name);
}

std::string SequenceAlignmentParameter::get_state() {
  return("n/a");
}

void SequenceAlignmentParameter::save_to_file(uint128_t gen, double l) {
  save_count += 1;
  msa->saveToFile(save_count, gen, l);
}

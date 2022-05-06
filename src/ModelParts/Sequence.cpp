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

static const state_element gap_indicator = -1;

// Sequence Alignment class.

SequenceAlignment::SequenceAlignment(std::string name, std::string msa_out, std::string subs_out, const States* states) : domain_name(name) {
  this->states = states->possible;
  n_states = states->n;

  state_element_encode = states->state_to_int;
  state_element_decode = states->int_to_state;

  seqs_out_file = msa_out;
  substitutions_out_file = subs_out;

  this->rare_threshold = env.get<double>("MCMC.rare_threshold");
}

void SequenceAlignment::add(std::string name, std::string sequence_str) {
  // Adds extant sequence to alignment. This will not be sampled during the MCMC.
  std::vector<state_element> enc = encode_sequence(sequence_str);
  taxa_names_to_sequences[name] = enc;
}

void SequenceAlignment::add(std::string name) {
  // Adds sequence to alignment, that WILL be sampled during MCMC.
  // This is for the ancestral nodes.
  std::vector<state_element> enc = {};
  taxa_names_to_sequences[name] = enc;
}

double** create_state_probability_vector(unsigned int n_cols, unsigned int n_states) {
  double** m = new double*[n_cols];
  for(unsigned int i = 0; i < n_cols; i++) {
    m[i] = new double[n_states];
    for(unsigned int j = 0; j < n_states; j++) {
      m[i][j] = 0.0;
    }
  }

  return(m);
}

void SequenceAlignment::add_base(std::string name, const IO::FreqSequence &seq) {
  add(name, sequenceAsStr_highestFreq(seq));

  prior_state_distribution[name] = create_state_probability_vector(seq.size(), n_states);

  int pos = 0;
  for(auto it = seq.begin(); it != seq.end(); ++it) {
    for(auto jt = it->begin(); jt != it->end(); ++jt) {
      if(jt->state != '-') {
	prior_state_distribution[name][pos][state_element_encode[std::string(1, jt->state)]] = jt->freq;
      }
    }
    pos++;
  }
}

void SequenceAlignment::print() {
  std::cout << "SEQUENCES" << std::endl;
  for(std::map<std::string, std::vector<state_element>>::iterator it = taxa_names_to_sequences.begin(); it != taxa_names_to_sequences.end(); ++it) {
    std::cout << ">" << it->first << "\n" << decode_state_element_sequence(it->second) << std::endl;
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
    buffer << ">" << it->first << "\n" << decode_state_element_sequence(it->second) << std::endl;
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
	subs_buffer << state_element_decode[anc] << pos << state_element_decode[dec] << " ";
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
  std::cout << "\tAttaching " << name << "[" << id << "] states to tree." << std::endl;
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
	  vector<state_element> dsNodeLseq = *(n->left->decendant->sequences[name]);
	  vector<state_element> dsNodeRseq = *(n->right->decendant->sequences[name]);
	  vector<state_element> p = find_parsimony(dsNodeLseq, dsNodeRseq);
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

    std::vector<state_element> sequence = taxa_names_to_sequences[(*node)->name];
    std::vector<bool> gaps(sequence.size(), false);

    marginal_state_distribution[(*node)->name] = create_state_probability_vector(n_columns, n_states);

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
std::vector<state_element> SequenceAlignment::encode_sequence(const std::string &sequence) {
  /*
   * Takes a string representation of a sequence and returns vector of integers.
   * Also tracks the gaps in the alignment.
   */
  std::vector<state_element> encoded_sequence(sequence.length());

  for (unsigned int site = 0; site < sequence.length(); site++) {
    std::string current_pos = sequence.substr(site, 1);
    try {
      encoded_sequence.at(site) = state_element_encode.at(current_pos);
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
int SequenceAlignment::n_cols() {
  return(n_columns);
}

std::string SequenceAlignment::decode_state_element(state_element c) {
  return(state_element_decode[c]);
}

std::string SequenceAlignment::decode_state_element_sequence(std::vector<state_element> &enc_seq) {
  std::string decoded_sequence;
  for(std::vector<state_element>::iterator it = enc_seq.begin(); it != enc_seq.end(); ++it) {
    decoded_sequence.append(decode_state_element(*it));
  }
  return(decoded_sequence);
}

std::vector<state_element> SequenceAlignment::find_parsimony(const std::vector<state_element> &s1, const std::vector<state_element> &s2) {
  std::vector<state_element> p = {};
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

// TODO positions should be passed as argument.
void SequenceAlignment::reset_to_base(std::string name, const std::list<unsigned int>& positions) {
  for(auto pos_it = positions.begin(); pos_it != positions.end(); pos_it++) {
    for(unsigned int j = 0; j < n_states; j++) {
      marginal_state_distribution[name][*pos_it][j] = prior_state_distribution[name][*pos_it][j];
    }
  } 
}

void SequenceAlignment::normalize_state_probs(TreeNode* node, unsigned int pos) {
  double normalize_total = 0.0;

  for(unsigned int i = 0; i < n_states; i++) {
    normalize_total += marginal_state_distribution[node->name][pos][i];
  }

  if(normalize_total != 0.0) {
    for(unsigned int i = 0; i < n_states; i++) {
      marginal_state_distribution[node->name][pos][i] /= normalize_total;
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

double SequenceAlignment::find_state_prob_given_dec_branch(BranchSegment* branch, state_element state_i, double* state_probs, std::vector<Valuable*> rv, double u, unsigned int pos) {
  /*
   * Find state probability given decendent branch
   * state_probs = the marginal posterior distribution of the state at the node below.
   */

  double prob = 0.0;
  double focal_domain_prob = 0.0;
  double alt_domain_prob = 1.0;

  float t_b = branch->distance;

  bool include_alt_domain = env.get<bool>("MCMC.include_alternative_domain");

  for(state_element state_j = 0; state_j < (state_element)n_states; state_j++) {
    double state_prob = state_probs[state_j];

    if(state_prob != 0.0) {
      // Likelihood contribution of all substitutions - including alternate domains.
      for(BranchSegment::iterator it = branch->begin(pos); it != branch->end(); it++) {
	std::string domain = (*it).first;

	if(domain == this->domain_name) {
	  double rate = rv[state_j]->get_value();
	  if(state_i != state_j) {
	    // Normal Substitution
	    focal_domain_prob = calc_substitution_prob(rate, t_b, u);
	  } else {
	    // No substitution - or possibly virtual.
	    focal_domain_prob = calc_no_substitution_prob(rate, t_b, u);
	  }
	} else if(include_alt_domain) {
	  // Subsitutions in non focal domain.
	  Substitution sub = (*it).second;
	  std::map<std::string, state_element> context = {{domain, sub.anc_state},
	  						  {this->domain_name, state_i}};

	  RateVector* rv = branch->get_hypothetical_rate_vector(domain, context, pos);

	  if(sub.occuredp and (sub.anc_state != sub.dec_state)) {
	    // Substitution including virtual substitutions.
	    alt_domain_prob *= calc_substitution_prob(rv->rates[sub.dec_state]->get_value(), t_b, u);
	  } else {
	   alt_domain_prob *= calc_no_substitution_prob(rv->rates[sub.anc_state]->get_value(), t_b, u);
	  }
	}	
      }
      prob += (state_prob * focal_domain_prob * alt_domain_prob);

      // Reset state-specific probability terms.
      focal_domain_prob = 0.0;
      alt_domain_prob = 1.0;
    }
  }

  return(prob);
}

double SequenceAlignment::find_state_prob_given_anc_branch(BranchSegment* branch, state_element state_j, double* state_probs, TreeNode* node, double u, unsigned int pos) {
  /*
   * Find state probability given ancestor branch.
   */

  double prob = 0.0;
  double focal_domain_prob = 0.0;
  double alt_domain_prob = 1.0;

  float t_b = branch->distance;

  bool include_alt_domain = env.get<bool>("MCMC.include_alternative_domain");

  for(state_element state_i = 0; state_i < (signed char)n_states; ++state_i) {
    double state_prob = state_probs[state_i];
    if(state_prob != 0.0) {
      for(BranchSegment::iterator it = branch->begin(pos); it != branch->end(); it++) {
	std::string domain = (*it).first;

	if(domain == this->domain_name) {
	  // Focal Domain
	  std::map<std::string, state_element> context = {{this->domain_name, state_i}};
	  RateVector* rv = node->up->get_hypothetical_rate_vector(domain_name, context, pos);

	  double rate = rv->rates[state_j]->get_value(); // i -> j rate.
	  if(state_i != state_j) {
	    // Normal Substitution.
	    focal_domain_prob = calc_substitution_prob(rate, t_b, u); // Probability of the substitution.
	  } else {
	    // No substition - possibly virtual.
	    focal_domain_prob = calc_no_substitution_prob(rate, t_b, u);
	  }

	} else if(include_alt_domain) {
	  // Alternative domains.
	  Substitution sub = (*it).second;
	  std::map<std::string, state_element> context = {{domain, sub.anc_state},
							  {this->domain_name, state_i}};
	  
	  RateVector* rv = branch->get_hypothetical_rate_vector(domain, context, pos);
	  if(sub.occuredp and (sub.anc_state != sub.dec_state)) {
	    //Substitution including virtual substitutions.
	    alt_domain_prob *= calc_substitution_prob(rv->rates[sub.dec_state]->get_value(), t_b, u);
	  } else {
	    alt_domain_prob *= calc_no_substitution_prob(rv->rates[sub.anc_state]->get_value(), t_b, u);
	  }
	}
      }
      prob += (state_prob * focal_domain_prob * alt_domain_prob);

      // Reset state-specific probability terms.
      focal_domain_prob = 0.0;
      alt_domain_prob = 1.0;
    }
  }

  return(prob);
}

void SequenceAlignment::find_marginal_at_pos(TreeNode* node, unsigned int pos, TreeNode* left_node, TreeNode* right_node, TreeNode* up_node) {
  double u = node->SM->get_u();
  double left_prob = 1.0;
  double right_prob = 1.0;
  double up_prob = 1.0;

  unsigned long extended_state;
  RateVector* rv;

  for(state_element state_i = 0; state_i < (signed char)n_states; state_i++) {
    // Most likely to be 0.0 so evaluated first.
    if(up_node != nullptr) {
      if(not taxa_names_to_gaps[up_node->name][pos]) {
	up_prob = find_state_prob_given_anc_branch(node->up, state_i, marginal_state_distribution[up_node->name][pos], node, u, pos);

	// Return if probability if 0.0.
	if(up_prob == 0.0) {
	  marginal_state_distribution[node->name][pos][state_i] = 0.0;
	  return;
	}
      }
    }

    // Contribution of left branch.
    if((left_node != nullptr) and not taxa_names_to_gaps[left_node->name].at(pos)) {
      std::map<std::string, state_element> context = {{domain_name, state_i}};
      rv = left_node->up->get_hypothetical_rate_vector(domain_name, context, pos);

      left_prob = find_state_prob_given_dec_branch(left_node->up, state_i, marginal_state_distribution[left_node->name][pos], rv->rates, u, pos);

      if(left_prob == 0.0) {
	marginal_state_distribution[node->name][pos][state_i] = 0.0;
	return;
      }
    }

    // Contribution of right branch.
    if((right_node != nullptr) and (not taxa_names_to_gaps[right_node->name][pos])) {
      std::map<std::string, state_element> context = {{domain_name, state_i}};
      rv = right_node->up->get_hypothetical_rate_vector(domain_name, context, pos);

      right_prob = find_state_prob_given_dec_branch(right_node->up, state_i, marginal_state_distribution[right_node->name][pos], rv->rates, u, pos);
      if(right_prob == 0.0) {
	marginal_state_distribution[node->name][pos][state_i] = 0.0;
	return;
      }
    }

    marginal_state_distribution[node->name][pos][state_i] = left_prob * right_prob * up_prob;
  }
}

void SequenceAlignment::find_state_probs_dec_only(TreeNode* node, std::list<unsigned int> positions) {
  /*
   * Finds the marginal posterior distribution for each position at a given node.
   * Only uses infomation from nodes below - used for upward recursion.
   * Assumes that node is not a tip.
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
	// Always a left node.
	find_marginal_at_pos(node, *pos, node->left->decendant, right_node, nullptr);
	normalize_state_probs(node, *pos);
      }
    }
  }
}

// TODO refactor.
void SequenceAlignment::find_state_probs_all(TreeNode* node, std::list<unsigned int> positions) {
  // NOTE assumes not a tip.
  std::string name = node->name;
  std::vector<bool> gaps = taxa_names_to_gaps[name];

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
      find_marginal_at_pos(node, *pos, node->left->decendant, right_node, up_node);
      normalize_state_probs(node, *pos);
    }
  }
}

// Second recursion.
void SequenceAlignment::update_state_probs(TreeNode* node, unsigned int pos, TreeNode* up_node) {
  // NOTE we can assume up_node is not a nullptr.
  double u = node->SM->get_u();
  double* state_probs = marginal_state_distribution[node->name][pos];

  for(state_element state_j = 0; state_j < (signed char)n_states; state_j++) {
    if(state_probs[state_j] != 0.0) {
      //std::cout << "update: " << (unsigned int)state_j << std::endl;
      state_probs[state_j] *= find_state_prob_given_anc_branch(node->up, state_j, marginal_state_distribution[up_node->name][pos], node, u, pos);
    }
  }
}

void SequenceAlignment::fast_update_state_probs_tips(TreeNode* node, unsigned int pos, TreeNode* up_node) {
  /*
   * Only valid for tips.
   * Equivilent of:
   * reset_to_base()
   * update_state_probs(node, *pos, node->up->ancestral);
   * NOTE we can assume there is an up node at a tip.
   */

  double u = node->SM->get_u();
  double* state_probs = marginal_state_distribution[node->name][pos];
  double* base_state_probs = prior_state_distribution[node->name][pos];

  for(state_element state_j = 0; state_j < (state_element)n_states; state_j++) {
    if(base_state_probs[state_j] != 0.0) {
      state_probs[state_j] = base_state_probs[state_j] * find_state_prob_given_anc_branch(node->up, state_j, marginal_state_distribution[up_node->name][pos], node, u, pos);
    } else {
      state_probs[state_j] = 0.0;
    }
  }
}

// Third Recursion
int SequenceAlignment::pick_state_from_probabilities(TreeNode* node, int pos) {
  /*
   * Picks a state from the marginal posterior distribution (taxa_names _to_state_probs).
   */
  double* probs = marginal_state_distribution[node->name][pos];

  //state_element e = 0;
  //std::cout << node->name << " " << pos << " ";

  //if(node->up != nullptr) {
  //  e = node->up->ancestral->sequences[domain_name]->at(pos);
  //  std::cout << (signed int)e << " ";
  //}

  double r = Random();
  double acc = 0.0;
  int val = -1;

  int largest_i = -1;
  double largest_val = 0.0;
  for(unsigned int i = 0; i < n_states; i++) {
    acc += probs[i];

    if(r < acc and val == -1) {

      if(probs[i] < this->rare_threshold) {
	// Skip rare events
	probs[i] = 0.0;
	continue;
      }

      val = i;
      probs[i] = 1.0;
    } else {
      if(probs[i] > largest_val) {
	largest_val = probs[i];
	largest_i = i;
      }
      probs[i] = 0.0;
    }
  }

  //if(val == -1) {
  // std::cerr << "Error: incorrectly picking a state: " <<  pos << " " << node->name << std::endl;
  // exit(EXIT_FAILURE);
  //}

  if(val == -1) {
    probs[largest_i] = 1.0;
    return(largest_i);
  } else {
    return(val);
  }
}

void SequenceAlignment::pick_states_for_node(TreeNode* node, const std::list<unsigned int>& positions) {
  std::vector<bool> gaps = taxa_names_to_gaps[node->name];

  for(auto pos = positions.begin(); pos != positions.end(); ++pos) {
    // Pick state from marginal distributions.
    if(gaps[*pos]) {
      taxa_names_to_sequences[node->name][*pos] = -1;
    } else {
      taxa_names_to_sequences[node->name][*pos] = pick_state_from_probabilities(node, *pos);
    }
  }
}

void SequenceAlignment::reconstruct_expand(const std::list<TreeNode*>& recursion_path, const std::list<unsigned int>& positions) {
  /*
   * This function needs a better name.
   * This recalculates the marginal posteriors of each node and picks sequences, which again alters the posterior
   * distribution.
   * Starts at a randomly chosen node and propagates outwards from there.
   */

  for(auto it = recursion_path.begin(); it != recursion_path.end(); it++) {
    TreeNode* node = *it;
    if(node->isTip()) {
      // Tip Node.
      //reset_to_base(node->name, positions);

      std::vector<bool> gaps = taxa_names_to_gaps[node->name];

      for(auto pos = positions.begin(); pos != positions.end(); ++pos) {
	if(not gaps[*pos]) {
	  //update_state_probs(node, *pos, node->up->ancestral);
	  fast_update_state_probs_tips(node, *pos, node->up->ancestral);
	  normalize_state_probs(node, *pos);
	}
      }
    } else {
      // Internal Node.
      find_state_probs_all(node, positions);
    }

    pick_states_for_node(node, positions);
  }
}

// SAMPLING AND RECURSION
void SequenceAlignment::reverse_recursion(const std::list<unsigned int>& positions) {
  /* 
   * Initial reverse recurstion to start marginal posterior calculations of states at each node.
   * Starts at tips and works up the tree to the root.
   * Nodes are ordered in the list such that they are visted in order up the tree.
   */ 

  const std::list<TreeNode*> nodes = this->tree->nodes();
  
  for(auto n = nodes.begin(); n != nodes.end(); ++n) {
    if(not (*n)->isTip()) {
      find_state_probs_dec_only(*n, positions);
    } else {
      // This is important as states at tips can be uncertain.
      reset_to_base((*n)->name, positions);
    }

    //std::cout << (*n)->name << " [ ";
    //for(signed int i = 0; i < (signed int)n_states; i++) {
    //  std::cout << marginal_state_distribution[(*n)->name][0][i] << " ";
    //} 
    //std::cout << "]" << std::endl;
  }
}

sample_status SequenceAlignment::sample_with_double_recursion(const std::list<unsigned int>& positions) {
  const std::list<TreeNode*> nodes = tree->nodes();

  reverse_recursion(positions);

  // 2nd Recursion - Reverse recursion.
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
      // Note does not call for root node.
      if((not gaps[*pos]) and (not (up_node == nullptr))) {
	update_state_probs(node, *pos, up_node);
	normalize_state_probs(node, *pos);
      }
    }

    pick_states_for_node(node, positions);
  }

  return(sample_status({false, true, true}));
}

sample_status SequenceAlignment::sample_with_triple_recursion(const std::list<unsigned int>& positions) {
  const std::list<TreeNode*> nodes = tree->nodes();

  reverse_recursion(positions);

  // 2nd Recursion - Reverse recursion.
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
      // Note does not call for root node.
      if((not gaps[*pos]) and (not (up_node == nullptr))) {
	update_state_probs(node, *pos, up_node);
	normalize_state_probs(node, *pos);
      }

    }
  }

  // 3rd Recursion - picking states.
  reconstruct_expand(tree->get_recursion_path(tree->rand_node()), positions);

  return(sample_status({false, true, true}));
}

// These functions are not critical they are useful though for other user who may not know how they are breaking simPLEX.
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

SequenceAlignmentParameter::SequenceAlignmentParameter(SequenceAlignment* msa, unsigned int n_sample) : SampleableComponent("SequenceAlignment") {
  save_count = -1;
  this->msa = msa;
  this->n_sample = n_sample;
  this->n_cols = msa->n_cols();

  // Options
  this->triple_recursion = env.get<bool>("MCMC.triple_recursion");

  // Exit program if invalid environment settings.
  if(n_sample < 1) {
    std::cerr << "Error: MCMC.position_sample_count must be greater than 0." << std::endl;
    exit(EXIT_FAILURE);
  }

  if(n_sample > n_cols) {
    std::cerr << "Error: cannot sample " << n_sample << " from alignment with " << n_cols << " columms." << std::endl;
    std::cerr << "Maximum value of MCMC.position_sample_count is " << n_cols << "." << std::endl;
    exit(EXIT_FAILURE);
  } else if(n_sample == n_cols) {
    this->sample_loc = 0;
  } else {
    this->sample_loc = rand() % n_cols;
  }
}

void SequenceAlignmentParameter::print() {
  std::cout << "SequenceAlignment-" << msa->domain_name << std::endl;
}

std::string SequenceAlignmentParameter::get_type() {
  return("SEQUENCE_ALIGNMENT");
}

sample_status SequenceAlignmentParameter::sample() {
  // Makes list of all positions.
  // Leaving room for a feature where not all positions are sampled each time.

  std::cout << "Sampling " << msa->domain_name << ": "<< sample_loc << "->";

  //Find the positions to be sampled.
  unsigned int last_pos = 0;
  std::list<unsigned int> positions = {};
  while(positions.size() < n_sample) {
    positions.push_back(sample_loc);
    last_pos = sample_loc;

    sample_loc++;
    if(sample_loc >= n_cols) {
      if(not (positions.size() == n_sample)) {
	std::cout << last_pos << ",0->";
      }
      sample_loc = 0;
    }
  }

  std::cout << last_pos << std::endl;

  if(this->triple_recursion) {
    // Triple recursion
    return(msa->sample_with_triple_recursion(positions));
    // Double recursion
  } else {
    return(msa->sample_with_double_recursion(positions));
  }
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

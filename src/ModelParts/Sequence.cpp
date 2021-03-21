#include <sstream>
#include <algorithm>

#include "Sequence.h"
#include "../Environment.h"
#include "../IO/Files.h"

#include "Trees/TreeParts.h"
#include "Trees/Tree.h"

extern double Random();
extern Environment env;
extern IO::Files files;

std::vector<std::string> aa({"A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"});
std::vector<std::string> nucleotides({"A", "T", "C", "G"});

// Sequence Alignment class.

SequenceAlignment::SequenceAlignment(std::string name, std::string msa_out, std::string subs_out, const States* states) : name(name) {
  this->states = states->possible;
  n_states = states->n;

  state_to_integer = states->state_to_int;
  integer_to_state = states->int_to_state;

  seqs_out_file = msa_out;
  substitutions_out_file = subs_out;
}

static const signed char gap_indicator = -1;

void SequenceAlignment::add(std::string name, std::string sequence_str) {
  // Adds extant sequence to alignment. This will not be sampled during the MCMC.
  std::vector<signed char> enc = EncodeSequence(sequence_str);
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
  add(name, sequenceAsStr(seq));

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
  seqs_out_identifier = name + "_sequences_out";
  files.add_file(seqs_out_identifier, seqs_out_file, IOtype::OUTPUT);

  substitutions_out_identifier = name + "_substitutions_out";
  files.add_file(substitutions_out_identifier, substitutions_out_file, IOtype::OUTPUT);
  files.write_to_file(substitutions_out_identifier, "I,GEN,LogL,Ancestral,Decendant,Substitutions\n");

  n_columns = (*taxa_names_to_sequences.begin()).second.size();
}

void SequenceAlignment::saveToFile(int save_count, int gen, double l) {
  std::ostringstream buffer;
  buffer << "#" << save_count << ":" << gen << ":" << l << std::endl;
  for(auto it = taxa_names_to_sequences.begin(); it != taxa_names_to_sequences.end(); ++it) {
    buffer << ">" << it->first << "\n" << decodeSequence(it->second) << std::endl;
  }

  files.write_to_file(seqs_out_identifier, buffer.str());

  std::ostringstream subs_buffer;

  std::list<BranchSegment*> branches = tree->get_branches();
  for(auto it = branches.begin(); it != branches.end(); ++it) {
    subs_buffer << save_count << "," << gen << "," << l << ",";
    subs_buffer << (*it)->ancestral->name << "," << (*it)->decendant->name << ",[ ";
    std::vector<Substitution> subs = (*it)->get_substitutions(name);
    for(unsigned int pos = 0; pos < subs.size(); pos++) {
      if(subs[pos].occuredp == true) {
	int anc = (*it)->ancestral->sequences[name]->at(pos);
	int dec = (*it)->decendant->sequences[name]->at(pos);
	if(anc != dec) {
	    // Virtual Substitution.
	    //subs_buffer << integer_to_state[anc] << pos << integer_to_state[dec] << "* ";
	  //} else {
	    // Normal Substitution.
	    subs_buffer << integer_to_state[anc] << pos << integer_to_state[dec] << " ";
	  }
      }
    }
    subs_buffer << "]\n";
  }
  files.write_to_file(substitutions_out_identifier, subs_buffer.str());
}

void SequenceAlignment::syncWithTree(std::string name, unsigned int id, Tree* tree) {
  this->tree = tree;
  std::cout << "Sync Hidden States: " << id << " " << name << std::endl;
  std::list<TreeNode*> nodes = tree->nodes();

  for(auto it = nodes.begin(); it != nodes.end(); ++it) {
    TreeNode* n = *it;

    //  if(id < n->hidden_state_sequences.size()) {
    //  std::cerr << "Error: syncing hidden states out of order with tree." << std::endl;
    //  exit(EXIT_FAILURE);
    //}

    //n->hidden_state_sequences.resize(id+1, nullptr);

    if(taxa_names_to_sequences.count(n->name)) {
      n->sequences[name] = &(taxa_names_to_sequences.at(n->name));
    } else {
      // NORMAL RUN
      // only tip sequences are needed.
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
      // This doesn't deal with branch nodes.
      //std::cout << "Calculate State Probabilites: " << name << std::endl;
      for(int pos = 0; pos < sequence.size(); pos++) {
	if(taxa_names_to_gaps[(*node)->left->decendant->name][pos]) {
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
std::vector<signed char> SequenceAlignment::EncodeSequence(const std::string &sequence) {
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

std::string SequenceAlignment::decodeChar(signed char &c) {
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

void SequenceAlignment::reset_base_probabilities(std::list<int> positions) {
  for(auto seq_name = base_sequences.begin(); seq_name != base_sequences.end(); ++seq_name) {
    //std::cout << *seq_name << " [ ";
    for(unsigned int i = 0; i < n_columns; i++) {
      // std::cout << "[ ";
      for(unsigned int j = 0; j < n_states; j++) {
	taxa_names_to_state_probs[*seq_name][i][j] = base_taxa_state_probs[*seq_name][i][j];
	//std::cout << taxa_names_to_state_probs[*seq_name][i][j] << " ";
      }
      //std::cout << "] ";
    }
    //std::cout << "]" << std::endl;
  }
}

void SequenceAlignment::calculate_state_probabilities_pos(TreeNode* node, unsigned int pos, TreeNode* left_node, TreeNode* right_node, TreeNode* up_node) {
  double u = node->SM->get_u();
  RateVector* rv;
  float left_prob = 0.0;
  float right_prob = 0.0;
  float up_prob = 0.0;

  double total = 0.0;
  for(int i = 0; i < n_states; i++) {
    // Contribution of left branch.
    if(left_node != nullptr) {
      std::vector<signed char> left_seq = taxa_names_to_sequences[left_node->name];
      std::vector<bool> left_gaps = taxa_names_to_gaps[left_node->name];

      if(not left_gaps[pos]) {
	float left_t_b = left_node->distance;
	
	left_prob = 0.0;
	unsigned long extended_state = left_node->get_hash_state_by_pos(pos);
	//extended_state[name] = i;
	rv = node->SM->selectRateVector({pos, i, name, extended_state});
	for(int j = 0; j < n_states; j++) {
	  double rate = rv->rates[left_seq.at(pos)]->get_value();
	  double state_prob = taxa_names_to_state_probs[left_node->name][pos][j];
	  left_prob += (state_prob * rate * left_t_b)/(1.0 + (u * left_t_b));
	}

	// Probability of staying the same.
	left_prob += taxa_names_to_state_probs[left_node->name][pos][i] / (1.0 + (u * left_t_b)); 
      } else {
	  left_prob = 1.0;
      }
    } else {
      left_prob = 1.0;
    }

    // Contribution of right branch.
    if(right_node != nullptr) {
      std::vector<signed char> right_seq = taxa_names_to_sequences[right_node->name];
      std::vector<bool> right_gaps = taxa_names_to_gaps[right_node->name];

      if(not right_gaps[pos]) {
	float right_t_b = right_node->distance;

	right_prob = 0.0;
	unsigned long extended_state = right_node->get_hash_state_by_pos(pos);
	//extended_state[name] = i; 
	rv = node->SM->selectRateVector({pos, i, name, extended_state});
	for(int j = 0; j < n_states; j++) {
	  double rate = rv->rates[right_seq.at(pos)]->get_value();
	  double state_prob = taxa_names_to_state_probs[right_node->name][pos][j];
	  right_prob += (state_prob * rate * right_t_b)/(1.0 + (u * right_t_b));
	}
	
	right_prob += taxa_names_to_state_probs[right_node->name][pos][i] / (1.0 + (u * right_t_b)); 
      } else {
	right_prob = 1.0;
      } 
    } else {
      right_prob = 1.0;
    }

    // Contribution of up branch.
    if(up_node != nullptr) {
      std::vector<signed char> up_seq = taxa_names_to_sequences[up_node->name];
      std::vector<bool> up_gaps = taxa_names_to_gaps[up_node->name];

      float up_t_b = node->distance;

      up_prob = 0.0;
      for(int j = 0; j < n_states; j++) {
	// Double check this.
	unsigned long extended_state = up_node->get_hash_state_by_pos(pos);
	//extended_state[name] = j;
	rv = node->SM->selectRateVector({pos, j, name, extended_state});
	double rate = rv->rates[i]->get_value();
	double state_prob = taxa_names_to_state_probs[up_node->name][pos][j];
	up_prob += (state_prob * rate * up_t_b)/(1.0 + (u * up_t_b));
      }

      // Probability of staying the same.
      up_prob += taxa_names_to_state_probs[up_node->name][pos][i] / (1.0 + (u * up_t_b)); 
    } else {
      up_prob = 1.0;
    }

    taxa_names_to_state_probs[node->name][pos][i] = left_prob * right_prob * up_prob;
    total += left_prob * right_prob * up_prob;

    if(isnan(left_prob * right_prob * up_prob)) {
      std::cerr << "Error: -nan" << std::endl;
      std::cerr << left_prob << " " << right_prob << " " << up_prob << std::endl;
      std::cerr << left_node->name << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  // Normalize
  if(total != 0.0) {
    for(int i = 0; i < n_states; i++) {
      taxa_names_to_state_probs[node->name][pos][i] /= total;
    }
  }
}

void SequenceAlignment::calculate_state_probabilities(TreeNode* node, std::list<int> positions) {
  std::string name = node->name;
  std::vector<bool> gaps = taxa_names_to_gaps[name];
  if(not node->isTip()) {
    // Left has always been sampled, as intermeidate nodes are connected by left and top.
    TreeNode* right_node;
    if(node->right) {
      right_node = node->right->decendant;
    } else {
      right_node = nullptr;
    }

    for(auto pos = positions.begin(); pos != positions.end(); ++pos) {
      if(not gaps[*pos]) {
	calculate_state_probabilities_pos(node, *pos, node->left->decendant, right_node, nullptr);
      }
    }
  }
}

int SequenceAlignment::pick_state_from_probabilities(TreeNode* node, int pos) {
  float* probs = taxa_names_to_state_probs[node->name][pos];

  //std::cout << node->name << ": [ ";
  //for(int i = 0; i < n_states; i++) {
  // std::cout << probs[i] << " ";
  //}
  //std::cout << "]" << std::endl;

  double r = Random();
  double acc = 0.0;
  int val = -1;
  for(int i = 0; i < n_states; i++) {
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

  std::list<int> positions = {};
  for(unsigned int i = 0; i < n_columns; i++) {
    positions.push_back(i);
  };

  for(auto n = nodes.begin(); n != nodes.end(); ++n) {
    (*n)->sampledp = false;
  }

  // Sample tip nodes if need be - this should be redundant now.
  reset_base_probabilities(positions);

  // Find state probabilities.
  for(auto n = nodes.begin(); n != nodes.end(); ++n) {
    if(not (*n)->isTip() and not (*n)->ready_to_sample()) {
      std::cerr << "Error: sampling of sequences is incorrectly ordered." << std::endl;
      exit(EXIT_FAILURE);
    } else {
      //std::cout << (*n)->name << " T: " << (*n)->isTip() << std::endl;
      (*n)->sampledp = true;
      calculate_state_probabilities(*n, positions);

      // if((*n)->isTip()) {
      //std::cout << "After: " << (*n)->name << " [ ";
      //for(int i = 0; i < n_columns; i++) {
      //  std::cout << "[ ";
      //  for(int j = 0; j < n_states; j++) {
      //    std::cout << taxa_names_to_state_probs[(*n)->name][i][j] << " ";
      //  }
      //  std::cout << "] ";
      //}
      //std::cout << "] " << std::endl;
      //}
    }
  }

  // Reverse recursion.
  for(auto n = nodes.rbegin(); n != nodes.rend(); ++n) {
    TreeNode* node = *n;
    std::vector<bool> gaps = taxa_names_to_gaps[node->name];

    // Reculaculate state probability vector - including up branch.
    // No need if root.
    if(true) {
      TreeNode* left_node;
      if(node->left) {
	left_node = node->left->decendant;
      } else {
	left_node = nullptr;
      }

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
	  calculate_state_probabilities_pos(node, *pos, left_node, right_node, up_node);
	}

	// If tip weight by data.
	if((*n)->isTip()) {
	  if(not gaps[*pos]) {
	    float total = 0.0;
	    for(int i = 0; i < n_states; i++) {
	      taxa_names_to_state_probs[node->name][*pos][i] *= base_taxa_state_probs[node->name][*pos][i];
	      total += taxa_names_to_state_probs[node->name][*pos][i];
	    }
	    
	    //Normalize.
	    for(int i = 0; i< n_states; i++) {
	      taxa_names_to_state_probs[node->name][*pos][i] /= total;
	    }
	  }
	}
      }
    }

    //Debug
    //if((*n)->isTip()) {
    // std::cout << "Final: " << (*n)->name << " [ ";
    // for(int i = 0; i < n_columns; i++) {
    //std::cout << "[ ";
    //for(int j = 0; j < n_states; j++) {
    //  std::cout << taxa_names_to_state_probs[(*n)->name][i][j] << " ";
    //}
    //std::cout << "] ";
    // }
    // std::cout << "] " << std::endl;
    //} 
  }

  // Pick state.
  for(auto n = nodes.rbegin(); n != nodes.rend(); ++n) {
    TreeNode* node = *n;
    std::vector<bool> gaps = taxa_names_to_gaps[node->name];
   
    for(auto pos = positions.begin(); pos != positions.end(); ++pos) {
      if(gaps[*pos]) {
	taxa_names_to_sequences[(*n)->name][*pos] = -1;
      } else {
	taxa_names_to_sequences[(*n)->name][*pos] = pick_state_from_probabilities(node, *pos);
      }
    }
  }

  // Reset all nodes such that the sampled flag is false.
  for(auto n = nodes.begin(); n != nodes.end(); ++n) {
    (*n)->sampledp = false;
  }

  return(sample_status({false, true, true}));
}

// This need to be much more robust - add better error handling.
bool SequenceAlignment::validate(std::list<std::string> seq_names, std::map<std::string, SequenceAlignment*> other_alignments) {
  for(auto n = seq_names.begin(); n != seq_names.end(); ++n) {
    if(taxa_names_to_sequences.find(*n) == taxa_names_to_sequences.end()) {
      std::cerr << "Error: sequence alignment " << name << " is missing sequence for " << *n << std::endl;
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
	std::cerr << "Error: in sequence alignment \"" << name << "\": sequences are not the same length as reference." << std::endl;
	exit(EXIT_FAILURE);
	return(false);
      } else {
	for(unsigned int i = 0; i < seq->second.size(); i++) {
	  if((seq->second[i] == -1) xor (it->second[i] == -1)) {
	    std::cout << seq->second[i] << " " << it->second[i] << std::endl;
	    std::cerr << "Error: in sequence alignment \"" << name << "\" in sequence " << seq->first << " at position " << i << " inconsistant gaps." << std::endl;
	    exit(EXIT_FAILURE);
	    return(false);
	  }
	} 
      }
    } else {
      std::cerr << "Error: in sequence alignment \"" << name << "\": sequence for \"" << seq->first << "\" is not found in reference." << std::endl;
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
  std::cout << "SequenceAlignment-" << msa->name << std::endl;
}

std::string SequenceAlignmentParameter::get_type() {
  return("SEQUENCE_ALIGNMENT");
}

sample_status SequenceAlignmentParameter::sample() {
  std::cout << "Sampling: " << msa->name << std::endl;
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

void SequenceAlignmentParameter::save_to_file(int gen, double l) {
  save_count += 1;
  msa->saveToFile(save_count, gen, l);
}

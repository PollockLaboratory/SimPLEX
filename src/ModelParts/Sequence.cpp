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

SequenceAlignment::SequenceAlignment(const States* states) {
  this->states = states->possible;
  n_states = states->n;

  state_to_integer = states->state_to_int;
  integer_to_state = states->int_to_state;
}

SequenceAlignment::SequenceAlignment(const SequenceAlignment &msa) {
  // The copy constructor.
  // Check whether this is truely copying.
  taxa_names_to_sequences = msa.taxa_names_to_sequences;
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

void SequenceAlignment::add_base(std::string name, std::string sequence_str) {
  base_sequences.push_front(name);
  add(name, sequence_str);
}

void SequenceAlignment::print() {
  std::cout << "SEQUENCES" << std::endl;
  for(std::map<std::string, std::vector<int>>::iterator it = taxa_names_to_sequences.begin(); it != taxa_names_to_sequences.end(); ++it) {
    std::cout << ">" << it->first << "\n" << decodeSequence(it->second) << std::endl;
  }
}

void SequenceAlignment::Initialize(std::list<SequenceAlignment*>* msa_List) {
  // Process sequences.
  MSA_list = msa_List;

  // Setup output.
  files.add_file("sequences_out", env.get<std::string>("OUTPUT.sequences_out_file"), IOtype::OUTPUT);

  n_columns = (*taxa_names_to_sequences.begin()).second.size();

}

void SequenceAlignment::Initialize(IO::RawMSA* &raw_msa) {
  for(auto it = raw_msa->seqs.begin(); it != raw_msa->seqs.end(); ++it) {
    add_base((*it).first, (*it).second);
  }

  // Setup output.
  files.add_file("sequences_out", env.get<std::string>("OUTPUT.sequences_out_file"), IOtype::OUTPUT);

  n_columns = (*taxa_names_to_sequences.begin()).second.size();
}

void SequenceAlignment::Initialize(IO::RawAdvMSA raw_msa) {
  for(auto it = raw_msa.seqs.begin(); it != raw_msa.seqs.end(); ++it) {
    add_base(it->first, sequenceAsStr(it->second));
  }

  // Add output file.

  n_columns = (*taxa_names_to_sequences.begin()).second.size();
}

void SequenceAlignment::saveToFile(int gen, double l) {
  static int i = -1;
  ++i;

  std::ostringstream buffer;
  buffer << "#" << i << ":" << gen << ":" << l << std::endl;
  for(auto it = taxa_names_to_sequences.begin(); it != taxa_names_to_sequences.end(); ++it) {
    buffer << ">" << it->first << "\n" << decodeSequence(it->second) << std::endl;
  }

  files.write_to_file("sequences_out", buffer.str());
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

void SequenceAlignment::syncWithTree(Tree* tree){
  this->tree = tree;

  std::list<TreeNode*> nodes = tree->nodes();
  for(auto it = nodes.begin(); it != nodes.end(); ++it) {
    TreeNode* n = *it;
    n->MSA = this;

    if(taxa_names_to_sequences.count(n->name)) {
      n->sequence = &(taxa_names_to_sequences.at(n->name));
    } else {
      if(env.ancestral_sequences) {
	// ANCESTRAL SEQUENCES KNOWN
	// In the case when ancestral sequences are known there can be no missing sequences.
	// New sequences should not be created.
	std::cerr << "Error: Missing sequence for \"" << n->name << "\"." << std::endl;
	exit(EXIT_FAILURE);
      } else {
	// NORMAL RUN
	// only tip sequences are needed.
	if(n->isTip()){
	  std::cerr << "Error: Missing sequence for \"" << n->name << "\"." << std::endl;
	  exit(EXIT_FAILURE);
	} else {
	  // Add new sequence to sequence alignments.
	  add(n->name);
	  n->sequence = &(taxa_names_to_sequences.at(n->name));
	  
	  // Fill missing sequences/
	  if(n->left != 0 and n->right == 0) {
	    // Internal Continous.
	    TreeNode* dsNode = n->left->decendant; // ds = downstream.
	    *(n->sequence) = *(dsNode->sequence);
	  } else {
	    // Root or internal branch.
	    vector<int> dsNodeLseq = *(n->left->decendant->sequence);
	    vector<int> dsNodeRseq = *(n->right->decendant->sequence);
	    vector<int> p = findParsimony(dsNodeLseq, dsNodeRseq);
	    *(n->sequence) = p;
	  }
	}
      }
    }
  }

  identify_gaps();
}

void SequenceAlignment::syncHiddenWithTree(unsigned int id, Tree* tree) {
  std::cout << "Sync Hidden States: " << id << std::endl;
  std::list<TreeNode*> nodes = tree->nodes();

  for(auto it = nodes.begin(); it != nodes.end(); ++it) {
    TreeNode* n = *it;

    if(id < n->hidden_state_sequences.size()) {
      std::cerr << "Error: syncing hidden states out of order with tree." << std::endl;
      exit(EXIT_FAILURE);
    }

    n->hidden_state_sequences.resize(id+1, nullptr);

    if(taxa_names_to_sequences.count(n->name)) {
      n->hidden_state_sequences[id] = &(taxa_names_to_sequences.at(n->name));
    } else {
      // NORMAL RUN
      // only tip sequences are needed.
      if(n->isTip()){
	std::cerr << "Error: Missing sequence for \"" << n->name << "\"." << std::endl;
	exit(EXIT_FAILURE);
      } else {
	// Add new sequence to sequence alignments.
	add(n->name);
	n->hidden_state_sequences[id] = &(taxa_names_to_sequences.at(n->name));
	  
	// Fill missing sequences/
	if(n->left != 0 and n->right == 0) {
	  // Internal Continous.
	  TreeNode* dsNode = n->left->decendant; // ds = downstream.
	  *(n->hidden_state_sequences[id]) = *(dsNode->hidden_state_sequences[id]);
	} else {
	  // Root or internal branch.
	  vector<int> dsNodeLseq = *(n->left->decendant->hidden_state_sequences[id]);
	  vector<int> dsNodeRseq = *(n->right->decendant->hidden_state_sequences[id]);
	  vector<int> p = findParsimony(dsNodeLseq, dsNodeRseq);
	  *(n->hidden_state_sequences[id]) = p;
	}
      }
    }
  }

  //identify_gaps();
}
					     
void SequenceAlignment::identify_gaps() {
  std::list<TreeNode*> nodes = tree->nodes();
  for(auto node = nodes.begin(); node != nodes.end(); ++node) {
    //std::cout << (*node)->name << std::endl;

    std::vector<int> sequence = taxa_names_to_sequences[(*node)->name];
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
std::vector<int> SequenceAlignment::EncodeSequence(const std::string &sequence) {
  /*
   * Takes a string representation of a sequence and returns vector of integers.
   * Also tracks the gaps in the alignment.
   */
  std::vector<int> encoded_sequence(sequence.length());

  for (unsigned int site = 0; site < sequence.length(); site++) {
    std::string current_pos = sequence.substr(site, 1);
    try {
      encoded_sequence.at(site) = state_to_integer.at(current_pos);
    } catch(const std::out_of_range& e) {
      std::cerr << "Error: state \"" << current_pos << "\" in sequence alignment is not recognised. " << std::endl;
      exit(EXIT_FAILURE);
    }
  }
  
  return(encoded_sequence);
}

// Utilities
int SequenceAlignment::numCols() {
  return(n_columns);
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

std::list<std::string> SequenceAlignment::getNodeNames() {
  std::list<std::string> names;
  for(auto it = taxa_names_to_sequences.begin(); it != taxa_names_to_sequences.end(); ++it) {
    names.push_back(it->first);
  }
  return(names);
}

// Ancestral Sequences - Ignore this.
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

void SequenceAlignment::sample_base_sequences(std::list<int> positions) {
  for(auto seq_name = base_sequences.begin(); seq_name != base_sequences.end(); ++seq_name) {
    for(auto pos = positions.begin(); pos != positions.end(); ++pos) {
      std::vector<int> seq = taxa_names_to_sequences[*seq_name];
      taxa_names_to_state_probs[*seq_name][*pos][seq.at(*pos)] = 1.0;
    }
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
      std::vector<int> left_seq = taxa_names_to_sequences[left_node->name];
      std::vector<bool> left_gaps = taxa_names_to_gaps[left_node->name];

      float left_t_b = left_node->distance;

      left_prob = 0.0;
      rv = node->SM->selectRateVector({pos, i});
      for(int j = 0; j < n_states; j++) {
	if(not left_gaps[pos]) {
	  double rate = rv->rates[left_seq.at(pos)]->get_value();
	  double state_prob = taxa_names_to_state_probs[left_node->name][pos][j];
	  left_prob += (state_prob * rate * left_t_b)/(1.0 + (u * left_t_b));
	} else {
	  left_prob = 1.0;
	}
      }

      // Probability of staying the same.
      left_prob += taxa_names_to_state_probs[left_node->name][pos][i] / (1.0 + (u * left_t_b)); 
    } else {
      left_prob = 1.0;
    }

    // Contribution of right branch.
    if(right_node != nullptr) {
      std::vector<int> right_seq = taxa_names_to_sequences[right_node->name];
      std::vector<bool> right_gaps = taxa_names_to_gaps[right_node->name];

      float right_t_b = right_node->distance;

      right_prob = 0.0;
      rv = node->SM->selectRateVector({pos, i});
      for(int j = 0; j < n_states; j++) {
	if(not right_gaps[pos]) {
	  double rate = rv->rates[right_seq.at(pos)]->get_value();
	  double state_prob = taxa_names_to_state_probs[right_node->name][pos][j];
	  right_prob += (state_prob * rate * right_t_b)/(1.0 + (u * right_t_b));
	} else {
	  right_prob = 1.0;
	}
      }

      right_prob += taxa_names_to_state_probs[right_node->name][pos][i] / (1.0 + (u * right_t_b)); 
    } else {
      right_prob = 1.0;
    }

    // Contribution of up branch.
    if(up_node != nullptr) {
      std::vector<int> up_seq = taxa_names_to_sequences[up_node->name];
      std::vector<bool> up_gaps = taxa_names_to_gaps[up_node->name];

      float up_t_b = node->distance;

      up_prob = 0.0;
      for(int j = 0; j < n_states; j++) {
	// Double check this.
	rv = node->SM->selectRateVector({pos, j});
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

  // Sample tip nodes if need be.
  sample_base_sequences(positions);

  // Find state probabilities.
  for(auto n = nodes.begin(); n != nodes.end(); ++n) {
    if(not (*n)->isTip() and not (*n)->ready_to_sample()) {
      std::cerr << "Error: sampling of sequences is incorrectly ordered." << std::endl;
      exit(EXIT_FAILURE);
    } else {
      (*n)->sampledp = true;
      calculate_state_probabilities(*n, positions); 
    }
  }

  // Reverse recursion.
  for(auto n = nodes.rbegin(); n != nodes.rend(); ++n) {
    TreeNode* node = *n;
    std::vector<bool> gaps = taxa_names_to_gaps[node->name];

    // Reculaculate state probability vector - including up branch.
    // No need if root.
    if(not (*n)->isTip()) {
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
	  calculate_state_probabilities_pos(node, *pos, node->left->decendant, right_node, up_node);
	}
      }

      // Don't need to go through nodes to update this.
      for(auto pos = positions.begin(); pos != positions.end(); ++pos) {
	if(gaps[*pos]) {
	  (*node->sequence)[*pos] = -1;
	} else {
	  (*node->sequence)[*pos] = pick_state_from_probabilities(node, *pos);
	}
      }
    }
  }

  // Reset all nodes such that the sampled flag is false.
  for(auto n = nodes.begin(); n != nodes.end(); ++n) {
    (*n)->sampledp = false;
  }

  return(sample_status({false, true, true}));
}

SequenceAlignmentParameter::SequenceAlignmentParameter(SequenceAlignment* msa) : SampleableComponent("SequenceAlignment") {
  this->msa = msa;
}

void SequenceAlignmentParameter::print() {
  std::cout << "SequenceAlignment" << std::endl;
}

std::string SequenceAlignmentParameter::get_type() {
  return("SEQUENCE_ALIGNMENT");
}

sample_status SequenceAlignmentParameter::sample() {
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
  msa->saveToFile(gen, l);
}

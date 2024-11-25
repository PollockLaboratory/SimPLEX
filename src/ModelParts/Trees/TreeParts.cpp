/*
 * TREE PARTS
 * This file holds the classes for the components of the tree:
 * 	- Tree nodes - Contains the sequences.
 * 	- Branche segments - Contains the branch lengths and location of substitutions.
 */

#include <iostream>
#include <math.h>

#include "TreeParts.h"
#include "../../Environment.h"

#include "../SubstitutionModels/RateVector.h"
#include "../SubstitutionModels/SubstitutionModel.h"

extern double Random();
extern Environment env;

// BRANCH SEGMENT
BranchSegment::BranchSegment(float distance) {
  this->distance = distance;
}

void BranchSegment::Initialize(unsigned int n_columns, std::list<std::string> state_domain_names) {
  this->n_pos = n_columns;

  for(auto it = state_domain_names.begin(); it != state_domain_names.end(); ++it) {
    rates[*it] = std::vector<RateVector*>(n_columns, NULL);
    substitutions[*it] = std::vector<Substitution>(n_columns, Substitution({false, 0, 0, nullptr}));
  }
}

BranchSegment::~BranchSegment() {
  decendant->up = 0;
}

std::ostream& operator<< (std::ostream &out, const BranchSegment &b) {
  out << b.distance;
  return out;
}

const Substitution& BranchSegment::get_substitution(std::string domain, unsigned int pos) {
  auto it = substitutions.find(domain);
  if(it == substitutions.end()) {
    std::cerr << "Error: domain not recognized in BranchSegment." << std::endl;
    exit(EXIT_FAILURE);
  }

  return(substitutions.at(domain)[pos]);
}

const std::vector<Substitution>& BranchSegment::get_substitutions(std::string domain) {
  auto it = substitutions.find(domain);
  if(it == substitutions.end()) {
    std::cerr << "Error: domain not recognized in BranchSegment." << std::endl;
    exit(EXIT_FAILURE);
  }
  return(substitutions.at(domain));
}

unsigned long BranchSegment::get_hypothetical_hash_state(std::map<std::string, state_element>& input_states, unsigned int pos) {
  /*
   * Finds the hypothetical hash_state for the context at the ancestral end of this branch.
   * Input states have priority over the state present at the ancestral sequences.
   * Uses substitutions to calculate hash_state.
   */
  
  std::map<std::string, state_element> context = {};
  for(auto it = substitutions.begin(); it != substitutions.end(); it++) {
    if(input_states.find(it->first) != input_states.end()) {
      state_element s = input_states[it->first];
      context[it->first] = s;
    } else {
      context[it->first] = substitutions[it->first][pos].anc_state;
    }
  }

  return(decendant->SM->get_hash_state(context));
}

RateVector* BranchSegment::get_hypothetical_rate_vector(std::string focal_domain, std::map<std::string, state_element>& context, unsigned int pos) {
  unsigned long extended_state = this->get_hypothetical_hash_state(context, pos);
  return(ancestral->SM->selectRateVector({focal_domain, extended_state}));
}

BranchSegment::iterator::iterator(BranchSegment& branch, unsigned int pos, bool end): branch(branch), pos(pos) {
  it = end ? branch.substitutions.end() : branch.substitutions.begin();
}

std::pair<std::string, Substitution> BranchSegment::iterator::operator*() const {
  return(std::pair<std::string, Substitution>(it->first, branch.get_substitution(it->first, pos)));
}

BranchSegment::iterator& BranchSegment::iterator::operator++(int) {
  it++;
  return(*this);
}

bool BranchSegment::iterator::operator!=(const iterator& rhs) const {
  return(this->it != rhs.it);
}

BranchSegment::iterator BranchSegment::begin(unsigned int pos) {
  return(BranchSegment::iterator(*this, pos, false));
}

BranchSegment::iterator BranchSegment::end() {
  return(BranchSegment::iterator(*this, 0, true));
}

// Key Statistics
inline void BranchSegment::update_rate_vectors() {
  std::vector<state_element>* seq = ancestral->sequences.begin()->second;
  for(unsigned int pos = 0; pos < seq->size(); pos++) {
    if((*seq)[pos] != -1) {
      unsigned long compound_state_hash = ancestral->get_hash_state(pos);

      // Set the rate vectors for each domain of the state.
      for(const auto& [state_domain, sequence] : ancestral->sequences) {
        if (this->ancestral->SM->is_static(state_domain)) continue;
        RateVector* rv = ancestral->SM->selectRateVector({ state_domain, compound_state_hash });
        if(rv == nullptr) {
          std::cerr << "Error: cannot find RateVector for position: " << pos << " " << state_domain << std::endl;
          exit(EXIT_FAILURE);
        }

        this->rates[state_domain][pos] = rv;
      }
    }
  }
}

void BranchSegment::set_new_substitutions() {
  // Set new substitutions for all state domains.
  for(auto& [state_domain, counts] : this->substitutions) {
    if (this->ancestral->SM->is_static(state_domain)) continue;

    std::vector<state_element> *anc_seq = (ancestral->sequences[state_domain]);
    std::vector<state_element> *dec_seq = (decendant->sequences[state_domain]);

    std::vector<RateVector*> rv_set = this->rates[state_domain];
    for(unsigned int pos = 0; pos < anc_seq->size(); pos++) {
      if(anc_seq->at(pos) == -1 or dec_seq->at(pos) == -1) {
        counts[pos] = { false, -1, -1, nullptr };
      } else {
        if(anc_seq->at(pos) != dec_seq->at(pos)) {
          // Normal substitutions.
          counts[pos] = {true, anc_seq->at(pos), dec_seq->at(pos), rv_set[pos]};
        } else {
          // No substitution - possibility of virtual substitution.
          // This could be faster - we just need the virtual substitution rate.
          double vir_rate = rv_set[pos]->rates[dec_seq->at(pos)]->get_value();
          double p = 1.0 - (1.0 / (1.0 + (vir_rate * distance)));
          //std::cout << vir_rate << " " << p << std::endl;
          if(Random() < p) {
            // Virtual Substitution.
            counts[pos] = {true, anc_seq->at(pos), dec_seq->at(pos), rv_set[pos]};
          } else {
            // No Virtual Substitution.
            counts[pos] = {false, anc_seq->at(pos), dec_seq->at(pos), rv_set[pos]};
          }
        }
      }
    }
  }
}

void BranchSegment::update() {
  this->update_rate_vectors();
  this->set_new_substitutions();
}

// TREE NODES

int TreeNode::unique_id = 0;

TreeNode::TreeNode(IO::RawTreeNode* raw_tree) {
  /* 
   * TreeNode Constructor - from rawTreeNode.
   */
  name = raw_tree->name;
  distance = raw_tree->distance;
  up = 0;
  left = 0;
  right = 0;
  sequences = {};
}

TreeNode::TreeNode(std::string n) {
  /* 
   * TreeNode Constructor.
   */
  name = n;
  up = 0;
  left = 0;
  right = 0;
  sequences = {};
}

TreeNode::TreeNode() {
  /* 
   * TreeNode Constructor - from name and distance.
   */
  name = "Node" + std::to_string(unique_id);
  unique_id++;
  up = 0;
  left = 0;
  right = 0;
  sequences = {};

  tagp = false;
}

TreeNode::~TreeNode() {
}

void TreeNode::connect_substitution_model(SubstitutionModel* sm) {
  SM = sm;
}

std::string TreeNode::toString() {
  std::string n = name + ":" + std::to_string(distance);
  if(left == 0 and right == 0) {
    return(n);
  } else if(left == 0) {
    return("(" + right->decendant->toString() + ")" + n);
  } else if(right == 0) {
    return("(" + left->decendant->toString() + ")" + n);
  } else {
    return("(" + left->decendant->toString() + "," + right->decendant->toString() + ")" + n);
  }
}

unsigned long TreeNode::get_hash_state(unsigned int pos) {
  /*
   * Returns the hash_state for a particular position in the sequences at this tree node.
   */
  std::map<std::string, state_element> context = {};

  // Adds all state domains.
  for(auto it = this->sequences.begin(); it != this->sequences.end(); ++it) {
    context[it->first] = (*(it->second))[pos];
  }

  return(SM->get_hash_state(context));
}

bool TreeNode::isTip() const {
  if(this->left == 0 and this->right == 0) {
    return(true);
  } else {
    return(false);
  }
}

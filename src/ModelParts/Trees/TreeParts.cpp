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

void BranchSegment::Initialize(unsigned int n_columns, std::map<std::string, std::list<std::string>> all_states) {
  rates = std::vector<RateVector*>(n_columns, NULL);
  substitutions = std::vector<Substitution>(n_columns, Substitution({false, 0, 0, nullptr}));

  for(auto it = all_states.begin(); it != all_states.end(); ++it) {
    hidden_rates[it->first] = std::vector<RateVector*>(n_columns, NULL);
    hidden_substitutions[it->first] = std::vector<Substitution>(n_columns, Substitution({false, 0, 0, nullptr}));
  }
}

BranchSegment::~BranchSegment() {
  decendant->up = 0;
}

std::ostream& operator<< (std::ostream &out, const BranchSegment &b) {
  out << b.distance;
  return out;
}

const std::vector<Substitution>& BranchSegment::get_substitutions() {
  return(substitutions);
}

// Key Statistics
inline void BranchSegment::update_rate_vectors() {
  std::vector<int>* seq = ancestral->sequence;
  for(unsigned int pos = 0; pos < seq->size(); pos++) {
    if((*seq)[pos] == -1) {
      rates[pos] = NULL;
    } else {
      std::map<std::string, int> states = ancestral->get_extended_state_by_pos(pos);
      rv_request rq = {pos, (*seq)[pos], "primary", states};
      rates[pos] = ancestral->SM->selectRateVector(rq);

      // Set hidden state rate vector - and primary.
      for(auto it = ancestral->hidden_state_sequences.begin(); it != ancestral->hidden_state_sequences.end(); ++it) {
	rv_request rq = {pos, (*(ancestral->hidden_state_sequences[it->first]))[pos], it->first, states};
	hidden_rates[it->first][pos] = ancestral->SM->selectRateVector(rq);
      }
    }
  }
}

void BranchSegment::set_new_substitutions() {
  double u = decendant->SM->get_u();

  std::vector<int> *anc = (ancestral->sequence);
  std::vector<int> *dec = (decendant->sequence);

  for(int pos = 0; pos < anc->size(); pos++) {
    // Primary.
    if(anc->at(pos) == -1 or dec->at(pos) == -1) {
      //substitutions[pos] = {false, anc[pos], dec[pos], rates[pos]};
    } else {
      if(anc->at(pos) != dec->at(pos)) {
	// Normal substitutions.
	substitutions[pos] = {true, anc->at(pos), dec->at(pos), rates[pos]};
      } else {
	// Possibility of virtual substitution.
	// Not sure this is exactly right.
	//float length = distance;

	double vir_rate = rates[pos]->rates[dec->at(pos)]->get_value();
	//double p = 1 - (1 / (1 + (vir_rate * distance)));
	double p = vir_rate / (1 - u + vir_rate);
	if(Random() < p) {
	  substitutions[pos] = {true, anc->at(pos), dec->at(pos), rates[pos]};
	} else {
	  substitutions[pos] = {false, anc->at(pos), dec->at(pos), rates[pos]};
	}
      }
    }

    // Hidden states here - and primary there is alot of redundancy now.
    for(auto domain = hidden_substitutions.begin(); domain != hidden_substitutions.end(); ++domain) {
      std::vector<int> *hidden_anc = (ancestral->hidden_state_sequences[domain->first]);
      std::vector<int> *hidden_dec = (decendant->hidden_state_sequences[domain->first]);
      
      for(unsigned int pos = 0; pos < hidden_anc->size(); pos++) {
	if(hidden_anc->at(pos) == -1 or hidden_dec->at(pos) == -1) {
	  hidden_substitutions[domain->first][pos] = {false, hidden_anc->at(pos), hidden_dec->at(pos), hidden_rates[domain->first][pos]};
	} else {
	  if(hidden_anc->at(pos) != hidden_dec->at(pos)) {
	    // Normal substitutions.
	    hidden_substitutions[domain->first][pos] = {true, hidden_anc->at(pos), hidden_dec->at(pos), hidden_rates[domain->first][pos]};
	  } else {
	    // Possibility of virtual substitution.
	    // Not sure this is exactly right.
	    //float length = distance;
	    
	    double vir_rate = hidden_rates[domain->first][pos]->rates[hidden_dec->at(pos)]->get_value();
	    //double p = 1 - (1 / (1 + (vir_rate * distance)));
	    double p = vir_rate / (1 - u + vir_rate);
	    if(Random() < p) {
	      hidden_substitutions[domain->first][pos] = {true, hidden_anc->at(pos), hidden_dec->at(pos), hidden_rates[domain->first][pos]};
	    } else {
	      hidden_substitutions[domain->first][pos] = {false, hidden_anc->at(pos), hidden_dec->at(pos), hidden_rates[domain->first][pos]};
	    }
	  }
	}
      }
    }
  }
}

void BranchSegment::update() {
  update_rate_vectors();
  set_new_substitutions();
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
  sampledp = false;
  hidden_state_sequences = {};
}

TreeNode::TreeNode(std::string n) {
  /* 
   * TreeNode Constructor.
   */
  name = n;
  up = 0;
  left = 0;
  right = 0;
  sampledp = false;
  hidden_state_sequences = {};
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
  sampledp = false;
  hidden_state_sequences = {};
}

TreeNode::~TreeNode() {
}

void TreeNode::connect_substitution_model(SubstitutionModel* sm) {
  SM = sm;
}

bool TreeNode::ready_to_sample() {
  bool rightp = true;
  if(right) {
    rightp = right->decendant->sampledp;
  }

  bool leftp = true;
  if(left) {
    leftp = left->decendant->sampledp;
  }
  return(leftp and rightp);
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

std::map<std::string, int> TreeNode::get_extended_state_by_pos(int pos) {
  std::map<std::string, int> states = {};
  states["primary"] = (*sequence)[pos];

  // Debug
  for(auto it = hidden_state_sequences.begin(); it != hidden_state_sequences.end(); ++it) {
    states[it->first] = (*(it->second))[pos];
  }

  return(states);
}

bool TreeNode::isTip() {
  if(this->left == 0 and this->right == 0) {
    return(true);
  } else {
    return(false);
  }
}

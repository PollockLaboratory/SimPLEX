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
  n_pos = n_columns;

  for(auto it = all_states.begin(); it != all_states.end(); ++it) {
    rates[it->first] = std::vector<RateVector*>(n_columns, NULL);
    substitutions[it->first] = std::vector<Substitution>(n_columns, Substitution({false, 0, 0, nullptr}));
  }
}

BranchSegment::~BranchSegment() {
  decendant->up = 0;
}

std::ostream& operator<< (std::ostream &out, const BranchSegment &b) {
  out << b.distance;
  return out;
}

const std::vector<Substitution>& BranchSegment::get_substitutions(std::string domain) {
  auto it = substitutions.find(domain);
  if(it == substitutions.end()) {
    std::cerr << "Error: domain not recognized in BranchSegment." << std::endl;
    exit(EXIT_FAILURE);
  }
  return(substitutions.at(domain));
}

// Key Statistics
inline void BranchSegment::update_rate_vectors() {
  std::vector<signed char>* seq = ancestral->sequences.begin()->second;
  for(unsigned int pos = 0; pos < seq->size(); pos++) {
    if((*seq)[pos] != -1) {
      unsigned long ex_state = ancestral->get_hash_state(pos);

      // Set the rate vectors for each of the states..
      for(auto it = ancestral->sequences.begin(); it != ancestral->sequences.end(); ++it) {
	rv_request rq = {pos, it->first, ex_state};
	
	if(ancestral->SM->selectRateVector(rq) == nullptr) {
	  std::cerr << "Error: cannot find RateVector for position: " << pos << " " << it->first << std::endl;
	  exit(EXIT_FAILURE);
	}

	rates[it->first][pos] = ancestral->SM->selectRateVector(rq);
      }
    }
  }
}

void BranchSegment::set_new_substitutions() {
  double u = decendant->SM->get_u();

  for(unsigned int pos = 0; pos < n_pos; pos++) {
    // Set new substitutions for all states..
    for(auto domain = substitutions.begin(); domain != substitutions.end(); ++domain) {
      std::vector<signed char> *anc_seq = (ancestral->sequences[domain->first]);
      std::vector<signed char> *dec_seq = (decendant->sequences[domain->first]);
      
      for(unsigned int pos = 0; pos < anc_seq->size(); pos++) {
	if(anc_seq->at(pos) == -1 or dec_seq->at(pos) == -1) {
	  substitutions[domain->first][pos] = {false, anc_seq->at(pos), dec_seq->at(pos), rates[domain->first][pos]};
	} else {
	  if(anc_seq->at(pos) != dec_seq->at(pos)) {
	    // Normal substitutions.
	    substitutions[domain->first][pos] = {true, anc_seq->at(pos), dec_seq->at(pos), rates[domain->first][pos]};
	  } else {
	    // Possibility of virtual substitution.
	    // Not sure this is exactly right.
	    //float length = distance;

	    double vir_rate = rates[domain->first][pos]->rates[dec_seq->at(pos)]->get_value();

	    //double p = 1 - (1 / (1 + (vir_rate * distance)));
	    double p = vir_rate / (1 - u + vir_rate);
	    if(Random() < p) {
	      substitutions[domain->first][pos] = {true, anc_seq->at(pos), dec_seq->at(pos), rates[domain->first][pos]};
	    } else {
	      substitutions[domain->first][pos] = {false, anc_seq->at(pos), dec_seq->at(pos), rates[domain->first][pos]};
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
  //std::map<std::string, int> states = {};

  // Adds all states.
  //for(auto it = sequences.begin(); it != sequences.end(); ++it) {
  // states[it->first] = (*(it->second))[pos];
  // }

  return(SM->get_hash_state(sequences, pos));
}

unsigned long TreeNode::get_hypothetical_hash_state(unsigned int pos, std::string domain, signed char state) {
  return(SM->get_hypothetical_hash_state(sequences, pos, domain, state));
}

bool TreeNode::isTip() {
  if(this->left == 0 and this->right == 0) {
    return(true);
  } else {
    return(false);
  }
}

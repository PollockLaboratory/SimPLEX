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
#include "../../SubstitutionCounts.h"

extern double Random();
extern Environment env;

// BRANCH SEGMENT

BranchSegment::BranchSegment(float distance) {
  this->distance = distance;
}

void BranchSegment::Initialize(unsigned int n_columns) {
  rates = std::vector<RateVector*>(n_columns, NULL);
  substitutions = std::vector<Substitution>(n_columns, Substitution({false, 0, 0, nullptr}));
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

double BranchSegment::get_rate(int pos, int dec_state) {
  double r = rates[pos]->rates[dec_state]->get_value();
  // Might not need this anymore.
  if(isnan(log(r))) {
    rates[pos]->rates[dec_state]->print();
  }
  return(r);
}

// Key Statistics
inline void BranchSegment::update_rate_vectors() {
  std::vector<int>* seq = ancestral->sequence;
  for(unsigned int pos = 0; pos < seq->size(); pos++) {
    if((*seq)[pos] == -1) {
      rates[pos] = NULL;
    } else {
      rv_request rq = {pos, (*seq)[pos]};
      rates[pos] = ancestral->SM->selectRateVector(rq);
    }
  }
}

void BranchSegment::set_new_substitutions() {
  double u = decendant->SM->get_u();

  std::vector<int> anc = *(ancestral->sequence);
  std::vector<int> dec = *(decendant->sequence);

  for(unsigned int pos = 0; pos < anc.size(); pos++) {
    if(anc[pos] == -1 or dec[pos] == -1) {
      substitutions[pos] = {false, anc[pos], dec[pos], rates[pos]};
    } else {
      if(anc[pos] != dec[pos]) {
	// Normal substitutions.
	substitutions[pos] = {true, anc[pos], dec[pos], rates[pos]};
      } else {
	// Possibility of virtual substitution.
	// Not sure this is exactly right.
	//float length = distance;
	double vir_rate = rates[pos]->rates[dec[pos]]->get_value();
	//double p = 1 - (1 / (1 + (vir_rate * distance)));
	double p = vir_rate / (1 - u + vir_rate);
	if(Random() < p) {
	  substitutions[pos] = {true, anc[pos], dec[pos], rates[pos]};
	} else {
	  substitutions[pos] = {false, anc[pos], dec[pos], rates[pos]};
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

std::string TreeNode::get_sequence() {
  std::string seq = "";
  for(unsigned int i = 0; i < sequence->size(); i++) {
    seq.append(MSA->integer_to_state[(*sequence)[i]]);
  }
  return(seq);
}

std::string TreeNode::state_at_pos(int i) {
  return(MSA->integer_to_state[(*sequence)[i]]);
}

bool TreeNode::isTip() {
  if(this->left == 0 and this->right == 0) {
    return(true);
  } else {
    return(false);
  }
}

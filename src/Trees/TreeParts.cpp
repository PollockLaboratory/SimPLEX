/*
 * TREE PARTS
 * This file holds the classes for the components of the tree:
 * 	- Tree nodes - Contains the sequences.
 * 	- Branche segments - Contains the branch lengths and location of substitutions.
 */

#include <iostream>

#include "TreeParts.h"
#include "../Environment.h"
#include "../SubstitutionCounts.h"

extern double Random();
extern Environment env;

// BRANCH SEGMENT

BranchSegment::BranchSegment(float distance) {
  this->distance = distance;
  rates = std::vector<RateVector*>(env.n, NULL);
  substitutions = std::vector<bool>(env.n, false);
}

BranchSegment::~BranchSegment() {
  decendant->up = 0;
}

std::ostream& operator<< (std::ostream &out, const BranchSegment &b) {
  out << b.distance;
  return out;
}

double BranchSegment::get_rate(int pos, int dec_state) {
  double r = rates[pos]->rates[dec_state]->getValue();
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
      rates[pos] = ancestral->SM->selectRateVector((*seq)[pos]);
    }
  }
}

void BranchSegment::set_new_substitutions() {
  double u = decendant->SM->get_u();
  for(auto it = substitutions.begin(); it != substitutions.end(); ++it) {
    if(Random() < u) {
      *it = true;
    } else {
      *it = false;
    }
  }
}

void BranchSegment::update() {
  update_rate_vectors();
}

void BranchSegment::update_counts(std::map<RateVector*, std::vector<int>>& subs_by_rateVector,
				  raw_counts& subs_by_branch) {
  std::vector<int> anc = *(ancestral->sequence);
  std::vector<int> dec = *(decendant->sequence);

  // Substitutions tracked.
  for(unsigned int pos = 0; pos < anc.size(); pos++) {
    if(substitutions[pos] == true and dec.at(pos) != -1) {
      // Adds both virtual substitutions and normal substitutions.
      subs_by_branch.num1subs += 1;
      subs_by_rateVector[rates[pos]][dec.at(pos)] += 1;
    } else {
      subs_by_branch.num0subs += 1;
    }
  } 
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
  sampled = false;
}

TreeNode::TreeNode(std::string n) {
  /* 
   * TreeNode Constructor.
   */
  name = n;
  up = 0;
  left = 0;
  right = 0;
  sampled = false;
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
  sampled = false;
}

// Sampling.
void TreeNode::sample_sequence() {
  if(right) {
    // Sample branch point node.
    for(unsigned int pos = 0; pos < sequence->size(); pos++) {
      if((*sequence)[pos] == -1) {
	continue;
      } else {
	if(left->substitutions[pos] == true or right->substitutions[pos]) {
	  (*sequence)[pos] = sample_single_position(pos);
	} else {
	  int left_dec = left->decendant->sequence->at(pos);
	  if(left_dec and right->decendant->sequence->at(pos)) {
	    (*sequence)[pos] = left_dec;
	  } else {
	    (*sequence)[pos] = sample_single_position(pos);
	  }
	}
      }
    }
  } else {
    // Only left decendant.
    for(unsigned int pos = 0; pos < sequence->size(); pos++) {
      // Skip sampling if gap.
      if((*sequence)[pos] == -1) {
	continue;
      } else {
	if(left->substitutions[pos] == true) {
	  (*sequence)[pos] = sample_single_position(pos);
	} else {
	  (*sequence)[pos] = left->decendant->sequence->at(pos);
	}
      }
    }
  }
}

TreeNode* TreeNode::sample() {
  // std::cout << "Name: " << name << " left: " << (bool)left << " right: " << (bool)right << std::endl;
  // Left has always been sampled, as intermediate nodes are connected by left and top.
  // Check right has been sampled.
  if(right) {
    if(not right->decendant->sampled) {
      return(this);
    }
  }

  sample_sequence();  
  sampled = true;

  if(up) {
    return(up->ancestral->sample());
  } else {
    // If at the root of the tree, return sampled node to end recurrsion.
    return(this);
  }
}

void branchLikelihood(double &l, int anc, int dec, float t_b, SubstitutionModel* SM) {
  /*
   * Calculates the likelihood of a branch, and adds it to the vector l.
   */
  // Doesn't take into account virtual substitutions.
  double u = SM->get_u();
  if(anc == dec) {
    l *= 1.0/(1.0 + t_b *u);
  } else {
    double rate = SM->selectRateVector(anc)->rates[dec]->getValue();
    l *= (rate*t_b)/(1.0 + t_b *u);
  }
}

std::vector<double> normalizeLikelihoods(std::vector<double> &l) {
  double total = 0.0;
  for(auto it = l.begin(); it != l.end(); ++it) {
    total += *it;
  }

  for(auto it = l.begin(); it != l.end(); ++it) {
    *it = *it / total;
  }
  return(l);
}

int TreeNode::sample_single_position(int pos) {
  int n = env.num_states;
  std::vector<double> l(n, 1.0);
  for(int state = 0; state < n; state++) {
    if(left) {
      int dec = left->decendant->sequence->at(pos);
      if(dec != -1) {
	branchLikelihood(l[state], state, dec, left->distance, SM);
      }
    }

    if(right) {
      int dec = right->decendant->sequence->at(pos);
      if(dec != -1) {
	branchLikelihood(l[state], state, dec, right->distance, SM);
      }
    }
  }

  l = normalizeLikelihoods(l);

  double r = Random();
  int i = 0;
  double c = l[0];

  while(r > c) {
    i++;
    c += l[i];
  }

  return(i);
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

bool TreeNode::isTip() {
  if(this->left == 0 and this->right == 0) {
    return(true);
  } else {
    return(false);
  }
}

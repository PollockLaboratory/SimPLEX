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
  rates = std::vector<RateVector*>(env.n, NULL);
  substitutions = std::vector<Substitution>(env.n, Substitution({false, 0, 0, nullptr}));
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
	float length = distance;
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

float** init_state_probabilities() {
  float** m = new float*[env.n];
  for(int i = 0; i < env.n; ++i) {
    m[i] = new float[env.num_states];
    for(int j = 0; j < env.num_states; j++) {
      m[i][j] = 0.0;
    }
  }
  return(m);
}

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

  state_probabilities = init_state_probabilities();
  gaps = new bool[env.n];
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

  state_probabilities = init_state_probabilities();
  gaps = new bool[env.n];
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

  state_probabilities = init_state_probabilities();
  gaps = new bool[env.n];
}

TreeNode* TreeNode::set_gaps() {
  if(isTip()) {
    sampledp = true;
    for(int i = 0; i < sequence->size(); i++) {
      if(sequence->at(i) == -1) {
	gaps[i] = true;
      } else {
	gaps[i] = false;
      }
    }
   
    return(up->ancestral);
  }

  if(not ready_to_sample()) {
    return(this);
  }

  // This doesn't deal with branch nodes.
  //std::cout << "Calculate State Probabilites: " << name << std::endl;
  for(int pos = 0; pos < env.n; pos++) {
    if(left->decendant->gaps[pos]) {
      if(right) {
	if(right->decendant->gaps[pos]) {
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

  sampledp = true;
 
  if(up) {
    // Returns pointer to node above to indicate it should be put on the queue ready to be sampled.
    return(up->ancestral);
  } else {
    // If at the root of the tree, return sampled node to end recursion.
    return(this);
  }
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

void TreeNode::calculate_state_probabilities_pos(int pos, TreeNode* left_node, TreeNode* right_node, TreeNode* up_node){
  double u = SM->get_u();
  float t_b = left->distance;
  RateVector* rv;
  float left_prob = 0.0;
  float right_prob = 0.0;
  float up_prob = 0.0;

  std::vector<int>* dec_seq = left_node->sequence;

  double total = 0.0;
  for(int i = 0; i < env.num_states; i++) {
    // Contribution of left branch.
    if(left_node != nullptr) {
      left_prob = 0.0;
      rv = SM->selectRateVector({pos, i});
      for(int j = 0; j < env.num_states; j++) {
	if(not left_node->gaps[pos]) {
	  double rate = rv->rates[dec_seq->at(pos)]->get_value();
	  left_prob += (left_node->state_probabilities[pos][j] * rate * t_b)/(1.0 + (u * t_b));
	} else {
	  left_prob = 1.0;
	}
      }

      // Probability of staying the same.
      left_prob += left_node->state_probabilities[pos][i] / (1.0 + (u * t_b)); 
    }

    // Contribution of right branch.
    if(right_node != nullptr) {
      //std::cout << "Right: " << right_node->name << std::endl;
      right_prob = 0.0;
      rv = SM->selectRateVector({pos, i});
      for(int j = 0; j < env.num_states; j++) {
	if(not right_node->gaps[pos]) {
	  double rate = rv->rates[right_node->sequence->at(pos)]->get_value();
	  right_prob += (right_node->state_probabilities[pos][j] * rate * t_b)/(1.0 + (u * t_b));
	} else {
	  right_prob = 1.0;
	}
      }

      // Probability of staying the same.
      if(isnan(right_node->state_probabilities[pos][i])) {
	std::cerr << "Error: another nan: " << right_node->name << std::endl;
	std::cerr << "[ ";
	for(int j = 0; j < env.num_states; j++) {
	  std::cerr << right_node->state_probabilities[pos][j] << " ";
	}
	std::cerr << "]" << std::endl;
	exit(EXIT_FAILURE);
      }
      
      right_prob += right_node->state_probabilities[pos][i] / (1.0 + (u * t_b)); 
    } else {
      right_prob = 1.0;
    }

    // Contribution of up branch.
    if(up_node != nullptr) {
      up_prob = 0.0;
      for(int j = 0; j < env.num_states; j++) {
	// Double check this.
	rv = SM->selectRateVector({pos, j});
	double rate = rv->rates[i]->get_value();
	up_prob += (up_node->state_probabilities[pos][j] * rate * t_b)/(1.0 + (u * t_b));
      }

      // Probability of staying the same.
      up_prob += up_node->state_probabilities[pos][i] / (1.0 + (u * t_b));
    } else {
      up_prob = 1.0;
    }
    
    state_probabilities[pos][i] = left_prob * right_prob * up_prob;
    total += left_prob * right_prob * up_prob;

    if(isnan(left_prob * right_prob * up_prob)) {
      std::cerr << "Error: -nan" << std::endl;
      std::cerr << left_prob << " " << right_prob << " " << up_prob << std::endl;
      std::cerr << left->decendant->name << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  // Normalize
  if(total != 0.0) {
    for(int i = 0; i < env.num_states; i++) {
      state_probabilities[pos][i] = state_probabilities[pos][i] / total;
    }
  }
}

TreeNode* TreeNode::calculate_state_probabilities(const std::list<int>& positions) {
  if(isTip()) {
    sampledp = true;
    for(auto pos = positions.begin(); pos != positions.end(); ++pos){
      if(not gaps[*pos]) {
	state_probabilities[*pos][sequence->at(*pos)] = 1.0;
      }
    }
 
    return(up->ancestral);
  }

  if(not ready_to_sample()) {
    return(this);
  }

  // Left has always been sampled, as intermediate nodes are connected by left and top.
  // Check right has been sampled.

  TreeNode* right_node;
  if(right) {
    right_node = right->decendant;
    if(not right->decendant->sampledp) {
      // Pointer to this node is returned indicating it should be put in the queue to be sampled later.
      return(this);
    } 
  } else {
    right_node = nullptr;
  }
  
  // This doesn't deal with branch nodes - not sure if this is true.
  //std::cout << "Calculate State Probabilites: " << name << std::endl;
  for(auto pos = positions.begin(); pos != positions.end(); ++pos){
    if(not gaps[*pos]) {
      calculate_state_probabilities_pos(*pos, left->decendant, right_node, nullptr);
    }
  }

  sampledp = true;
 
  if(up) {
    // Returns pointer to node above to indicate it should be put on the queue ready to be sampled.
    return(up->ancestral);
  } else {
    // If at the root of the tree, return sampled node to end recursion.
    return(this);
  }
}

// Picking the state from the pre-calculated probabilites.

int TreeNode::pick_state_from_probabilities(int pos) {
  float* probs = state_probabilities[pos];

  double r = Random();
  double acc = 0.0;
  int val = -1;
  for(int i = 0; i < env.num_states; i++) {
    acc += probs[i];
    if(r < acc and val == -1) {
      val = i;
      probs[i] = 1.0;
    } else {
      probs[i] = 0.0;
    }
  }

  if(val == -1) {
    std::cerr << "Error: incorrectly picking a state: " <<  pos << " " << name << std::endl;
    std::cerr << "[ ";
    for(int pos = 0; pos < env.n; pos++) {
      std::cerr << "[" << pos << "- " << gaps[pos] << " " << right->decendant->gaps[pos] << " " << left->decendant->gaps[pos] << " ]";
    }
    std::cerr << "]" << std::endl;
    exit(EXIT_FAILURE);
  }

  return(val);
}

void TreeNode::pick_sequences(const std::list<int>& positions) {
  if(isTip()) {
    return;
  }

  TreeNode* right_node;
  if(right) {
    right_node = right->decendant;
  } else {
    right_node = nullptr;
  }

  // Recalculate state probability vector - including the up branch.
  // No need if at the root.
  if(up != nullptr) {
    for(auto pos = positions.begin(); pos != positions.end(); ++pos){
      if(not gaps[*pos]) {
	  calculate_state_probabilities_pos(*pos, left->decendant, right_node, up->ancestral);
      }
    }
  }

  for(auto pos = positions.begin(); pos != positions.end(); ++pos){
    if(gaps[*pos]) {
      (*sequence)[*pos] = -1;
    } else {
      (*sequence)[*pos] = pick_state_from_probabilities(*pos);
    }
  }

  // Recursively call the remainder of the tree.
  if(right) {
    right->decendant->pick_sequences(positions);
  }

  if(left) {
    left->decendant->pick_sequences(positions);
  }
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

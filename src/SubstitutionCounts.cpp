#include "SubstitutionCounts.h"

#include "Environment.h"
#include "IO/Files.h"
#include "ModelParts/Trees/Tree.h"

extern Environment env;
extern IO::Files files;

SubstitutionCounts::SubstitutionCounts() {
}

SubstitutionCounts::SubstitutionCounts(std::vector<RateVector*> rvs, std::list<float> b_lens) {
  subs_by_rateVector = {};
  subs_by_branch = {};
  for(auto it = rvs.begin(); it != rvs.end(); ++it) {
    subs_by_rateVector[*it] = std::vector<int>(env.num_states, 0);
  }

  for(auto it = b_lens.begin(); it != b_lens.end(); ++it) {
    subs_by_branch[*it] = {0,0};
  }
}

void SubstitutionCounts::print() {
  // By Rate Vector.
  std::cout << "Substitutions by Rate Vector:" << std::endl;
  for(auto it = subs_by_rateVector.begin(); it != subs_by_rateVector.end(); ++it) {
    std::cout << "[\t";
    for(auto jt = it->second.begin(); jt != it->second.end(); ++jt) {
      std::cout << *jt << "\t";
    }
    std::cout << "] - " << it->first->name << std::endl;
  }

  // By Branch length.
  std::cout << "Substitutions by Branch Length:" << std::endl;
  for(auto it = subs_by_branch.begin(); it != subs_by_branch.end(); ++it) {
    std::cout << "[ 0:" << it->second.num0subs << "\t1: " << it->second.num1subs << "\t] - " << it->first << std::endl;
  }
}

CountsParameter::CountsParameter(SubstitutionCounts* counts) : AbstractComponent("SubstitutionCounts."), counts(counts) {
}

void CountsParameter::link_to_tree(TreeParameter* tp) {
  this->tree = tp->get_tree_ptr();
  this->add_dependancy(tp);
}

void CountsParameter::refresh() {
  *counts = SubstitutionCounts(tree->SM->get_RateVectors(), tree->get_branch_lengths());
  tree->update_counts(*counts);
}

void CountsParameter::print() {
  std::cout << "CountsParameter." << std::endl;
}

double CountsParameter::record_state(int gen, double l) {
  return(0.0);
}

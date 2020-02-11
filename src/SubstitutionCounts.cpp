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
    std::cout << "] - " << it->first->get_name() << std::endl;
  }

  // By Branch length.
  std::cout << "Substitutions by Branch Length:" << std::endl;
  for(auto it = subs_by_branch.begin(); it != subs_by_branch.end(); ++it) {
    std::cout << "[ 0:" << it->second.num0subs << "\t1: " << it->second.num1subs << "\t] - " << it->first << std::endl;
  }
}

// CountsParameter

std::ofstream CountsParameter::out_file;

CountsParameter::CountsParameter(SubstitutionCounts* counts, TreeParameter* tp) : AbstractComponent("SubstitutionCounts."), counts(counts) {

  this->tree = tp->get_tree_ptr();
  this->add_dependancy(tp);

  files.add_file("substitution_counts", env.get<std::string>("OUTPUT.counts_out_file"), IOtype::OUTPUT);
  out_file = files.get_ofstream("substitution_counts");

  // Print csv header to outfile.
  out_file << "RateVector,State";
  RateVector* rv = tree->get_SM()->get_RateVectors().front();
  for(int i = 0; i < rv->size(); i++) {
    out_file << "," << rv->get_state_by_pos(i);
  }
  out_file << std::endl;
}

void CountsParameter::refresh() {
  *counts = SubstitutionCounts(tree->get_SM()->get_RateVectors(), tree->get_branch_lengths());

  const std::list<BranchSegment*> branchList = tree->get_branches();
  for(auto it = branchList.begin(); it != branchList.end(); ++it) {
    BranchSegment* b = *it;
    for(auto sub = b->get_substitutions().begin(); sub != b->get_substitutions().end(); ++sub) {
      if(sub->occuredp == true) {
	// Adds both virtual substitutions and normal substitutions.
	counts->subs_by_branch[b->distance].num1subs += 1;
	counts->subs_by_rateVector[sub->rate_vector][sub->dec_state] += 1;
      } else {
	counts->subs_by_branch[b->distance].num0subs += 1;
      }
    }
  }
}

void CountsParameter::print() {
  std::cout << "CountsParameter." << std::endl;
}

double CountsParameter::record_state(int gen, double l) {
  for(auto it = counts->subs_by_rateVector.begin(); it != counts->subs_by_rateVector.end(); ++it) {
    out_file << it->first->get_name() << "," << it->first->get_state();
    for(auto jt = it->second.begin(); jt != it->second.end(); ++jt) {
      out_file << "," << *jt;
    }
    out_file << std::endl;
  }
  return(0.0);
}

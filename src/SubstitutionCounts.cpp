#include "SubstitutionCounts.h"

#include "Environment.h"
#include "IO/Files.h"

#include "ModelParts/Trees/Tree.h"
#include "ModelParts/SubstitutionModels/RateVector.h"

#include <sstream>
#include <string>

extern Environment env;
extern IO::Files files;

SubstitutionCounts::SubstitutionCounts() {
}

SubstitutionCounts::SubstitutionCounts(std::vector<RateVector*> rvs, std::list<float> b_lens, std::map<std::string, States> all_states) {
  // Make empty structures ready for counts.
  subs_by_rateVector = {};
  subs_by_branch = {};

  for(auto it = rvs.begin(); it != rvs.end(); ++it) {
    subs_by_rateVector[*it] = std::vector<int>((*it)->rates.size(), 0);
  }

  for(auto it = b_lens.begin(); it != b_lens.end(); ++it) {
    subs_by_branch[*it] = {0,0};
  }
}

void SubstitutionCounts::clear() {
  // Clear structures ready for counts.
  for(auto it = subs_by_rateVector.begin(); it != subs_by_rateVector.end(); ++it) {
    for(auto jt = it->second.begin(); jt != it->second.end(); ++jt) {
     *jt = 0;
    }
  }

  for(auto it = subs_by_branch.begin(); it != subs_by_branch.end(); it++) {
    it->second.num0subs = 0;
    it->second.num1subs = 0;
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

//std::ofstream CountsParameter::out_file;

CountsParameter::CountsParameter(SubstitutionCounts* counts, Tree* tree, std::map<std::string, std::list<std::string>> all_states) : AbstractComponent("SubstitutionCounts."), counts(counts), all_states(all_states) {
  hidden = true;

  //this->tree = tp->get_tree_ptr();
  //this->add_dependancy(tp);
  this->tree = tree;

  files.add_file("counts_by_ratevector_out", env.get<std::string>("OUTPUT.counts_out_file"), IOtype::OUTPUT);

  // Print csv header to outfile.
  std::ostringstream buffer;
  buffer << "I,GEN,LogL,RateVector,State";
  RateVector* rv = tree->get_SM()->get_RateVectors().front();
  for(int i = 0; i < rv->size(); i++) {
    buffer << "," << rv->get_state_by_pos(i);
  }
  buffer << ",Total" << std::endl;

  files.write_to_file("counts_by_ratevector_out", buffer.str());
}

void CountsParameter::fix() {
}

void CountsParameter::refresh() {
  std::cout << "Refreshing counts" << std::endl;
  // Create new structs for counts.
  static bool updated_table = false;
  if(not updated_table) {
    *counts = SubstitutionCounts(tree->get_SM()->get_RateVectors(), tree->get_branch_lengths(), tree->SM->get_all_states());
    updated_table = true;
  }

  counts->clear();

  std::map<std::string, States> all_state_domains = tree->SM->get_all_states();

  // Track counts by branch segment length and rate vector.
  const std::list<BranchSegment*> branchList = tree->get_branches();

  for(auto it = branchList.begin(); it != branchList.end(); ++it) {
    BranchSegment* b = *it;
    for(auto jt = all_state_domains.begin(); jt != all_state_domains.end(); ++jt) {
      std::string domain = jt->first;
      for(auto sub = b->get_substitutions(domain).begin(); sub != b->get_substitutions(domain).end(); ++sub) {
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
}

void CountsParameter::print() {
  std::cout << "CountsParameter." << std::endl;
}

std::string CountsParameter::get_state_header() {
  return(name);
}

std::string CountsParameter::get_state() {
  return("n/a");
}

void CountsParameter::save_to_file(boost::multiprecision::uint128_t gen, double l) {
  static int index = -1;
  index++;

  std::ostringstream buffer;

  for(auto it = counts->subs_by_rateVector.begin(); it != counts->subs_by_rateVector.end(); ++it) {
    buffer << index << "," << gen << "," << l << "," << it->first->get_name() << "," << it->first->get_state();
    unsigned int total = 0;
    for(auto jt = it->second.begin(); jt != it->second.end(); ++jt) {
      total += *jt;
      buffer << "," << *jt;
    }
    buffer << "," << total << std::endl;
  }

  files.write_to_file("counts_by_ratevector_out", buffer.str());
}

std::string CountsParameter::get_type() {
  return("COUNTS_PARAMETER");
}

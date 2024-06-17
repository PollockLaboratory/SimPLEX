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

SubstitutionCounts::SubstitutionCounts(std::vector<RateVector*> rvs, std::list<float> b_lens) {
  this->base_virtual = env.get<int>("MCMC.base_virtual");
  // Make empty structures ready for counts.
  subs_by_rateVector = {};
  subs_by_branch = {};

  for(auto it = rvs.begin(); it != rvs.end(); ++it) {
    subs_by_rateVector[*it] = std::vector<double>((*it)->rates.size(), 0.0);
  }

  for(auto it = b_lens.begin(); it != b_lens.end(); ++it) {
    subs_by_branch[*it] = {0,0};
  }
}

void SubstitutionCounts::clear() {
  // Clear structures ready for counts.
  for(auto it = subs_by_rateVector.begin(); it != subs_by_rateVector.end(); ++it) {
    for(unsigned int i = 0; i < it->second.size(); i++ ) {
      if((int)i == (it->first)->state) {
	// Always 1 virtual subtitution.
	(it->second)[i] = this->base_virtual;
      } else {
	(it->second)[i] = 0.0;
      }
    }
  }

  for(auto it = subs_by_branch.begin(); it != subs_by_branch.end(); it++) {
    it->second.num0subs = 0.0;
    it->second.num1subs = 0.0;
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
CountsParameter::CountsParameter(SubstitutionCounts* counts, Tree* tree, std::map<std::string, std::list<std::string>> all_states) : AbstractComponent("SubstitutionCounts."), counts(counts), all_states(all_states) {
  hidden = true;

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
  // Create new structs for counts.
  static bool updated_table = false;
  if(not updated_table) {
    *counts = SubstitutionCounts(tree->get_SM()->get_RateVectors(), tree->get_branch_lengths());
    updated_table = true;
  }

  this->counts->clear();

  std::map<std::string, States> all_state_domains = tree->SM->get_all_states();

  // Track counts by branch segment length and rate vector.
  const std::list<BranchSegment*> branchList = tree->get_branches();

  // DEBUG
  std::map<std::string, std::vector<int>> sub_counts = {};
  std::map<std::string, double> vir_sub_counts = {};

  // RESET COUNTS
  for(const auto& [state_domain, _] : all_state_domains) {
    int n_cols = (*branchList.begin())->get_substitutions(state_domain).size();
    sub_counts[state_domain] = std::vector<int>(n_cols, 0);
    vir_sub_counts[state_domain] = 0.0;
  }

  //for(auto it = branchList.begin(); it != branchList.end(); ++it) {
  for(const auto& branch : tree->get_branches()) {
    for(const auto& [state_domain, _] : all_state_domains) {
      if (tree->get_SM()->is_static(state_domain)) continue;
      //for(auto jt = all_state_domains.begin(); jt != all_state_domains.end(); ++jt) {
      int pos = 0;
      for(auto sub = branch->get_substitutions(state_domain).begin(); sub != branch->get_substitutions(state_domain).end(); ++sub) {
        if(sub->dec_state == -1) {
          // Skip gaps.
          pos++;
          continue;
        }

        if(sub->anc_state != sub->dec_state) {
          // Normal substitutions.
          counts->subs_by_branch[branch->distance].num1subs += 1;
          counts->subs_by_rateVector[sub->rate_vector][sub->dec_state] += 1;
          sub_counts[state_domain][pos] += 1;
          //std::cout << state_domain << " " << (unsigned int)sub->anc_state << " " << (unsigned int)sub->dec_state
          //          << " " << branch->ancestral->name
          //          << " " << branch->decendant->name
          //          << " " << pos
          //          << std::endl;
        } else {
          // Virtual substitutions.
          // Adds the expected virtual substitution count. Unlikely to be integer.
          double virtual_subs = sub->rate_vector->rates[sub->anc_state]->get_value() * branch->distance;
          vir_sub_counts[state_domain] += virtual_subs;

          counts->subs_by_branch[branch->distance].num1subs += virtual_subs;
          counts->subs_by_branch[branch->distance].num0subs += (1.0 - virtual_subs);
          counts->subs_by_rateVector[sub->rate_vector][sub->dec_state] += virtual_subs;
        }
        
        pos++;
      }
    }
  }

  //DEBUG
  //for(auto it = sub_counts.begin(); it != sub_counts.end(); ++it) {
  //  int sum = 0;
  //  std::cout << it->first << " [ ";
  //  for(auto jt = it->second.begin(); jt != it->second.end(); ++jt) {
  //    sum += *jt;
  //    std::cout << *jt << " ";
  //  }
  //  std::cout << "] " << sum << " " << vir_sub_counts[it->first]
	//      << " " << (float)vir_sub_counts[it->first]/(sum+vir_sub_counts[it->first]) << std::endl;
  //}
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

  for (const auto& [rate_vector, pcounts] : counts->subs_by_rateVector) {
    buffer << index << "," << gen << "," << l << "," << rate_vector->get_name() << "," << rate_vector->get_state();
    unsigned int total = 0;
    for(auto jt = pcounts.begin(); jt != pcounts.end(); ++jt) {
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

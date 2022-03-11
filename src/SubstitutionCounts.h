#ifndef SubstitutionCounts_h_
#define SubstitutionCounts_h_

#include <map>
#include <vector>
#include <boost/multiprecision/cpp_int.hpp>

#include "ModelParts/AbstractComponent.h"
#include "ModelParts/SubstitutionModels/States.h"

class RateVector;

struct raw_counts {
  int num0subs = 0;
  int num1subs = 0;
};

class SubstitutionCounts {
 public:
  SubstitutionCounts();
  SubstitutionCounts(std::vector<RateVector*>, std::list<float>, std::map<std::string, States>);

  std::map<RateVector*, std::vector<int>> subs_by_rateVector;
  std::map<float, raw_counts> subs_by_branch;
  void print();
};

class Tree;
class AncestralStatesParameter;

class CountsParameter : public AbstractComponent {
private:
  static std::ofstream out_file;
  SubstitutionCounts* counts;
  Tree* tree;
  std::map<std::string, std::list<std::string>> all_states;
public:
  CountsParameter(SubstitutionCounts*, Tree*, std::map<std::string, std::list<std::string>>);
  void fix() override;
  void refresh() override;
  void print() override;

  std::string get_state_header() override;
  std::string get_state() override;

  std::string get_type() override;

  void save_to_file(boost::multiprecision::uint128_t gen, double l);
};

#endif

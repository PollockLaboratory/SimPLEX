#ifndef SubstitutionCounts_h_
#define SubstitutionCounts_h_

#include "ModelParts/SubstitutionModels/RateVector.h"
#include "ModelParts/AbstractComponent.h"

struct raw_counts {
  int num0subs = 0;
  int num1subs = 0;
};

class SubstitutionCounts {
 public:
  SubstitutionCounts();
  SubstitutionCounts(std::vector<RateVector*>, std::list<float>);

  std::map<RateVector*, std::vector<int>> subs_by_rateVector;
  std::map<float, raw_counts> subs_by_branch;
  void print();
};

class Tree;
class TreeParameter;

class CountsParameter : public AbstractComponent {
private:
  static std::ofstream out_file;
  SubstitutionCounts* counts;
  Tree* tree;
public:
  CountsParameter(SubstitutionCounts*, TreeParameter*);
  void refresh() override;
  void print() override;
  double record_state(int gen, double l) override;
  std::string get_type() override;
};

#endif

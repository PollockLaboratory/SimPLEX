#ifndef SubstitutionCounts_h_
#define SubstitutionCounts_h_

#include "SubstitutionModels/Components/RateVector.h"

class SubstitutionCounts {
 public:
  SubstitutionCounts();
  SubstitutionCounts(std::vector<RateVector*>, std::list<float>);

  std::map<RateVector*, std::vector<int>> subs_by_rateVector;
  std::map<float, std::pair<int, int>> subs_by_branch; // First = num0subs, Second = num1subs
  void print();
};

#endif

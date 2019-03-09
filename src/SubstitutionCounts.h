#ifndef SubstitutionCounts_h_
#define SubstitutionCounts_h_

#include "SubstitutionModels/Components/RateVector.h"

struct raw_counts {
  int num0subs = 0;
  int num1subs = 0;
};

// Minor

class SubstitutionCounts {
 public:
  SubstitutionCounts();
  SubstitutionCounts(std::vector<RateVector*>, std::list<float>);

  std::map<RateVector*, std::vector<int>> subs_by_rateVector;
  std::map<float, raw_counts> subs_by_branch; // First = num0subs, Second = num1subs
  void print();
};

#endif

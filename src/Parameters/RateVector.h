#ifndef RateVector_h_
#define RateVector_h_

#include <vector>
#include <set>
#include <unordered_set>
#include <functional>
#include <string>
#include <list>

#include "AbstractValue.h"
#include "ContinuousFloat.h"
#include "VirtualSubstitutionRate.h"

class BranchSegment; // Defined in Trees/Types/TreeParts.h

struct bpos {
  // Represents the position and branch that a rate vector applies to.
  BranchSegment* branch;
  int pos;

  bool operator<(const bpos& x) const {
    if(branch == x.branch) {
      if(pos < x.pos) {
	return(true);
      } else {
	return(false);
      }
    }
    if(branch < x.branch) return(true);
  }

  bool operator==(const bpos& x) const {
    return(pos == x.pos and branch == x.branch);
  }
};

namespace std {
  template<>
  struct hash<bpos> {
    inline size_t operator()(const bpos& x) const {
      std::hash<int> int_hash;
      std::hash<BranchSegment*> bs_hash;
      return(int_hash(x.pos) + bs_hash(x.branch));
    }
  };
}


// Fundamental collection type of parameters in substitution model.
class RateVector {
 public:
  RateVector(std::string, int state, std::vector<AbstractValue*>);

  std::string name;
  std::vector<AbstractValue*> rates;
  int state; // Determines the (ancestral) state that this rate vector applies to.

  void remove_location(int pos, BranchSegment* bs);
  void add_location(int pos, BranchSegment* bs);
  std::unordered_set<bpos> get_locations();
  void clear_locations();

  void update_counts();
  double get_logLikelihood();

  void print();
  std::unordered_set<bpos> locations; // List of branches that rate vector applies to.
 private:
  std::vector<int> counts;
  int size;
  static int IDc;
};

// Collections of rate vectors.
class RateVectorSet {
  // This collection holds all of the rate vectors currently in the substitution model.
 public:
  RateVectorSet();
  std::vector<RateVector*> c;
  void Initialize();

  RateVector*& operator[] (const int i);

  void add(RateVector* v);
  void get_counts();
  void clear_locations();
  void check_duplicate_locations();

  void print();
  void saveToFile(int gen, double l);
 private:
  static std::ofstream out_file;
};
#endif

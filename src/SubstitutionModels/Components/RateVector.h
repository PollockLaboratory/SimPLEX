#ifndef RateVector_h_
#define RateVector_h_

#include <vector>
#include <set>
#include <unordered_set>
#include <functional>
#include <string>
#include <list>
#include <map>

#include "AbstractComponent.h"
#include "AbstractValueTypes.h"
#include "SampleableValueTypes.h"

class BranchSegment; // Defined in Trees/Types/TreeParts.h

// Fundamental collection type of parameters in substitution model.
class RateVector {
 public:
  RateVector(std::string, int state, std::vector<AbstractValue*>);

  std::string name;
  std::vector<AbstractValue*> rates;
  int state; // Determines the (ancestral) state that this rate vector applies to.
  float operator[](int);
  float get_rate_ratio(int i);

  void print();

  std::map<int, int> valueID_to_state; // Maps a values ID to the state it applies to in the rates vector.
 private:
  void update_counts();
  std::vector<int> counts; // The counts of substitutions that apply to each rate.
  std::vector<double> logLikelihoods; // The logLikelihoods associated with each rate. count * log(rate).
  int size;
  static int IDc;
};

// Collections of rate vectors.
class RateVectorSet {
  // This collection holds all of the rate vectors currently in the substitution model.
 public:
  RateVectorSet();
  std::vector<RateVector*> col;
  void Initialize();

  RateVector*& operator[] (const int i);

  void add(RateVector* v);

  void print();
  void saveToFile(int gen, double l);
 private:
  static std::ofstream out_file;
};

#endif

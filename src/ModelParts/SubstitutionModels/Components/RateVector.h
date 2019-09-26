#ifndef RateVector_h_
#define RateVector_h_

#include <vector>
#include <set>
#include <unordered_set>
#include <functional>
#include <string>
#include <list>
#include <map>

#include "../../AbstractComponent.h"
#include "AbstractValueTypes.h"
#include "SampleableValueTypes.h"
#include "../../../IO/SubstitutionModelParser.h"

struct rv_request {
  // Struct representing a request for a rate vector.
  // More specific requests later, for example branch position.
  int pos;
  int state;
};

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
  const int& getID();

  void print();

  std::map<int, int> valueID_to_state; // Maps a values ID to the state it applies to in the rates vector.
 private:
  int id;
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

  void add(RateVector* v, IO::rv_use_class);
  void organize(int seqLen, int numStates);

  RateVector* select(rv_request);

  void print();
  void saveToFile(int gen, double l);
 private:
  std::map<int, IO::rv_use_class> id_to_uc;
  // Tree structure of Rate Vectors.
  // Organized via positions -> states -> possible RateVectors.
  std::vector<std::vector<std::list<RateVector*>>> rv_tree; 
  static std::ofstream out_file;
};

#endif

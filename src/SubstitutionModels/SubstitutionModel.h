#ifndef SubstitutionModel_h_
#define SubstitutionModel_h_

#include <fstream>
#include <sstream>

#include <list>
#include <set>
#include <vector>
#include <map>
#include <string>
#include <queue>

#include "Components/AbstractComponent.h"
#include "Components/ComponentSet.h"
#include "Components/RateVector.h"
#include "../IO/SubstitutionModelParser.h"

// Temparily defined in Component set before it finds a better home.
//struct States {
// int n;
// std::set<std::string> possible;
// std::map<std::string, int> state_to_int;
// std::map<int, std::string> int_to_state;
//};

class SubstitutionModel {
  class iterator;
 public:
  SubstitutionModel();
  void from_raw_model(IO::raw_substitution_model*);

  // States.
  States states;
  void add_state(std::string);
  void print_states();
  const States* get_states();

  // Rate Vectors.
  void organizeRateVectors(int seqLen, int numStates);
  RateVector* selectRateVector(rv_request);
  std::vector<RateVector*> get_RateVectors();
  SubstitutionModel::iterator changed_vectors_begin();

  // Parameters.
  UniformizationConstant* u;
  const double& get_u();

  void printParameters();
  int getNumberOfParameters();
  void get_current_parameters(std::list<rv_loc>&); // I think this is redundant now.

  // Sample.
  bool SampleParameters();
  void accept(); //After a model is sampled it must be accepted or rejected before next sampling.
  void reject();

  void saveToFile(int gen, double l);
  void Terminate();
 private:
  std::ofstream* substitution_model_out;
  void finalize();

  ComponentSet components;
  RateVectorSet rateVectors;

  class iterator:public std::iterator<std::output_iterator_tag, std::pair<RateVector*, int>> {
  public:
    explicit iterator(SubstitutionModel&, bool);
    const rv_loc& operator*() const;
    iterator& operator++();
    bool operator!=(const iterator &) const;
    bool at_end() const;
  private:
    inline void step_to_next_location();
    inline bool step_to_next_component();
    SubstitutionModel& sub_model;
    bool endQ;
    std::list<AbstractComponent*> changed_comps; // List of the components that have changes with recent sampling.
    std::queue<AbstractComponent*> cq;
    std::list<rv_loc>::iterator location; // The location within a rate vector that has changed.
    std::list<rv_loc>::iterator location_iter_end;
  };
};

#endif

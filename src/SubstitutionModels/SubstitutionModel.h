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

struct States {
  int n;
  std::set<std::string> possible;
  std::map<std::string, int> state_to_int;
  std::map<int, std::string> int_to_state;
};

class SubstitutionModel {
  class iterator;
 public:
  SubstitutionModel();
  virtual void Initialize() = 0;
  virtual void Initialize(int number_of_sites, std::vector<std::string> states) = 0;

  // States.
  States states;
  void add_state(std::string);
  void print_states();
  const States* get_states();

  // Rate Vectors.
  RateVector* selectRateVector(int state);
  void add_rate_vector(RateVector* v);
  std::vector<RateVector*> get_RateVectors();
  SubstitutionModel::iterator changed_vectors_begin();

  // Parameters.
  void printParameters();
  int getNumberOfParameters();
  void get_current_parameters(std::list<std::pair<RateVector*, int>>&);

  // Sample.
  bool SampleParameters();
  void accept(); //After a model is sampled it must be accepted or rejected before next sampling.
  void reject();

  void saveToFile(int gen, double l);
  virtual void Terminate();
 protected:
  void finalize();
 private:
  std::ofstream* substitution_model_out;

  ComponentSet components;
  RateVectorSet rateVectors;

  class iterator:public std::iterator<std::output_iterator_tag, std::pair<RateVector*, int>> {
   public:
    explicit iterator(SubstitutionModel&, bool);
    const std::pair<RateVector*, int>& operator*() const;
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
    std::list<std::pair<RateVector*, int>>::iterator location; // The location within a rate vector that has changed.
    std::list<std::pair<RateVector*, int>>::iterator location_iter_end;
  };
};

#endif

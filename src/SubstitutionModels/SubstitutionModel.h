#ifndef SubstitutionModel_h_
#define SubstitutionModel_h_

#include <fstream>
#include <sstream>

#include <list>
#include <vector>
#include <map>
#include <string>
#include <queue>

#include "Components/AbstractComponent.h"
#include "Components/ComponentSet.h"
#include "Components/RateVector.h"

class SubstitutionModel {
  class iterator;
 public:
  SubstitutionModel();
  virtual void Initialize(int number_of_sites, std::vector<std::string> states) = 0;

  RateVector* selectRateVector(int state);

  bool SampleParameters();
  void accept(); //After a model is sampled it must be accepted or rejected before next sampling.
  void reject();

  void printParameters();
  int getNumberOfParameters();
  void get_current_parameters(std::list<std::pair<RateVector*, int>>&);

  void get_counts();

  std::vector<RateVector*> get_RateVectors();

  void saveToFile(int gen, double l);
  virtual void Terminate();

  SubstitutionModel::iterator changed_vectors_begin();
 protected:
  void add_rate_vector(RateVector* v);
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
    iterator& operator++(int);
    bool operator!=(const iterator &) const;
    bool at_end() const;
  private:
    inline void step_to_next_location();
    inline bool step_to_next_component();
    bool endQ;
    SubstitutionModel& sub_model;
    std::list<AbstractComponent*> changed_comps; // List of the components that have changes with recent sampling.
    std::queue<AbstractComponent*> cq;
    std::list<std::pair<RateVector*, int>>::iterator location; // The location within a rate vector that has changed.
    std::list<std::pair<RateVector*, int>>::iterator location_iter_end;
  };
};

#endif

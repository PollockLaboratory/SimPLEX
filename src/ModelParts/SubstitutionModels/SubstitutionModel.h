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

#include "../AbstractComponent.h"
#include "../ComponentSet.h"
#include "RateVector.h"
#include "../../IO/SubstitutionModelParser.h"
#include "States.h"

class SubstitutionModel {
  class iterator;
public:
  SubstitutionModel(Valuable* u);

  // Reading IO.
  RateVector* create_rate_vector(IO::raw_rate_vector rv, Valuable* u);
  void from_raw_model(IO::raw_substitution_model*);

  // States.
  const States* get_states();

  // Rate Vectors.
  void organizeRateVectors(int seqLen, int numStates);
  RateVector* selectRateVector(rv_request);
  std::vector<RateVector*> get_RateVectors();
  SubstitutionModel::iterator modified_begin(AbstractComponent*);

  // Parameters.
  Valuable* u;
  const double& get_u();

  std::list<AbstractComponent*> get_all_parameters();

  void saveToFile(int gen, double l);
private:
  void configure_States(std::set<std::string>);
  void configure_RateVectors(std::list<IO::raw_rate_vector>);
  void configure_HiddenStates(std::map<std::string, std::set<std::string>>);

  States states;
  std::map<std::string, States> hidden_states;

  RateVectorSet rateVectors;

  // Iterator
  // Given a AbstractComponent will return an iterator to the start of all the rate vector locations that are modified.
  class iterator:public std::iterator<std::output_iterator_tag, std::pair<RateVector*, int>> {
  public:
    explicit iterator(SubstitutionModel&, bool, AbstractComponent*);
    const rv_loc& operator*() const;
    iterator& operator++();
    bool operator!=(const iterator &) const;
    bool at_end() const;
  private:
    inline void step_to_next_location();
    inline bool step_to_next_component();
    SubstitutionModel& sub_model;
    bool endQ;
    //std::list<SampleableValue*> changed_comps; // List of the components that have changes with recent sampling.
    std::list<rv_loc>::const_iterator location; // The location within a rate vector that has changed.
    std::list<rv_loc>::const_iterator location_iter_end;
    std::list<Valuable*>::const_iterator current_parameter;
    std::list<Valuable*>::const_iterator valuables_end;
  };
};

#endif

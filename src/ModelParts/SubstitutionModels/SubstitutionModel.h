#ifndef SubstitutionModel_h_
#define SubstitutionModel_h_

#include <list>
#include <set>
#include <map>
#include <string>

#include "../AbstractComponent.h"
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
  const States* get_states(std::string domain);
  std::map<std::string, States> get_all_states();

  // Rate Vectors.
  void organizeRateVectors();
  RateVector* selectRateVector(rv_request);
  std::vector<RateVector*> get_RateVectors();
  unsigned long get_hash_state(const std::map<std::string, std::vector<signed char>*>& sequences, unsigned int pos);
  unsigned long get_hypothetical_hash_state(const std::map<std::string, std::vector<signed char>*>& sequences, unsigned int pos, std::string domain_name, signed char state);
  
  SubstitutionModel::iterator modified_begin(AbstractComponent*);

  // Parameters.
  Valuable* u;
  const double& get_u();

  std::list<AbstractComponent*> get_all_parameters();
  
  void saveToFile(int gen, double l, std::map<RateVector*, std::vector<int>> counts_by_rv);
private:
  void configure_States(std::map<std::string, std::list<std::string>>);
  void configure_RateVectors(std::list<IO::raw_rate_vector>);

  RateVectorSet rateVectors;
  std::map<std::string, States> all_states;

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

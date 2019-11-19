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

struct States {
  int n; // Careful about indels.
  std::set<std::string> possible;
  std::map<std::string, int> state_to_int;
  std::map<int, std::string> int_to_state;
};

class SubstitutionModel {
  class iterator;
public:
  SubstitutionModel(Valuable* u);

  // Reading IO.
  SampleableValue* realize_component(IO::raw_param*);
  SampleableValue* retreive_component(int id);
  RateVector* create_rate_vector(States states, IO::raw_rate_vector rv, Valuable* u);
  void from_raw_model(IO::raw_substitution_model*);

  std::map<int, SampleableValue*> realized_params;
  std::map<int, IO::raw_param*> raw_params;
    
  // States.
  States states;
  void add_state(std::string);
  void print_states();
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
  void finalize();

  RateVectorSet rateVectors;

  // Iterator
  // Given a AbstractComponent will return an iterator to all the rate vector locations that are modified.
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
    std::list<SampleableValue*> changed_comps; // List of the components that have changes with recent sampling.
    std::queue<AbstractComponent*> cq;
    std::list<rv_loc>::iterator location; // The location within a rate vector that has changed.
    std::list<rv_loc>::iterator location_iter_end;
  };
};

#endif

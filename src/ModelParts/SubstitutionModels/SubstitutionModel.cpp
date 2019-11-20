#include <iostream>

#include "SubstitutionModel.h"
#include "../../Environment.h"
#include "../../IO/Files.h"

extern Environment env;
extern IO::Files files;

SubstitutionModel::SubstitutionModel(Valuable* u) : u(u) {
  states.n = 0;
  realized_params = {};
}

void SubstitutionModel::add_state(std::string s) {
  // Check if state is gap.
  states.possible.insert(s);
  states.state_to_int[s] = states.n;
  states.int_to_state[states.n] = s;
  states.n++;
}

SampleableValue* SubstitutionModel::realize_component(IO::raw_param* param) {
  SampleableValue* p;
  if(realized_params.find(param->ID) == realized_params.end()) {
    if(param->t == IO::FLOAT) {
      IO::raw_ContinuousFloat* float_param = dynamic_cast<IO::raw_ContinuousFloat*>(param);
      p = new ContinuousFloat(float_param->name, float_param->init, float_param->step_size, 0.0);
    } else if(param->t == IO::CATEGORY) {
      IO::raw_DiscreteFloat* cat_param = dynamic_cast<IO::raw_DiscreteFloat*>(param);
      // Need to sort duplicate rate categories out here.
      RateCategories* cats = new RateCategories("cats", cat_param->categories);
      p = new DiscreteFloat(cat_param->name, cats);
    } else {
      std::cerr << "Error: unknown parameter type for parameter \"" << param->name << "\"." << std::endl;
      exit(EXIT_FAILURE);
    }
  } else {
    p = realized_params[param->ID];
  }
  return(p);
}

SampleableValue* SubstitutionModel::retreive_component(int id) {
  if(realized_params.find(id) == realized_params.end()) {
    // Note: doing this on one line adds element to map before function call.
    SampleableValue* s = realize_component(raw_params[id]);
    realized_params[id] = s;
  }

  return(realized_params[id]);
}

RateVector* SubstitutionModel::create_rate_vector(States states, IO::raw_rate_vector rv, Valuable* u) {
  std::vector<Valuable*> rates(states.n, nullptr);
  int s = states.state_to_int[rv.uc.state];
  VirtualSubstitutionRate* vir_rate = new VirtualSubstitutionRate("tmp_name_virtual_substitution.", u); 
  for(int i = 0; i < states.n; i++) {
    IO::raw_param* param = rv.rates.front();
    if(param == nullptr) {
      std::cerr << "Error: nullptr for parameter in rate vector " << rv.name << " at position " << i << "." << std::endl;
      exit(EXIT_FAILURE);
    }
    if(i != s) {
      rates[i] = retreive_component(param->ID);
      // Add dependancies.
      vir_rate->add_rate(rates[i]);
    } else {
      // Virtual Substitution rate.
      vir_rate->name = param->name;
      rates[i] = vir_rate;
    }
    rv.rates.pop_front();
  }

  if(not rv.rates.empty()) {
    std::cerr << "Error: unexpected number of rates - " << states.n << " is expected (given the number of states), however " << rv.rates.size() << " additional rate is found." << std::endl;
    exit(EXIT_FAILURE);
  }

  RateVector* new_rv = new RateVector(rv.name, states.state_to_int[rv.uc.state], rates);
  return(new_rv);
}

void SubstitutionModel::from_raw_model(IO::raw_substitution_model* raw_sm) {
  raw_params = raw_sm->get_map_parameters();

  // Read in states.
  for(auto it = raw_sm->states.begin(); it != raw_sm->states.end(); ++it) {
	add_state(*it);
  }

  // Read in rate vectors.
  for(auto raw_rv = raw_sm->rv_list.begin(); raw_rv != raw_sm->rv_list.end(); ++raw_rv) {
    RateVector* rv = create_rate_vector(states, *raw_rv, u);
    rateVectors.add(rv, (*raw_rv).uc);
  }

  finalize();
}

void SubstitutionModel::print_states() {
 std::cout << "States: ";
 for(auto it = states.possible.begin(); it != states.possible.end(); ++it) {
   std::cout << *it << ":" << states.state_to_int[*it] << " ";
 }
 std::cout << "- n = " << states.possible.size() << std::endl;
}

const States* SubstitutionModel::get_states() {
  return(&states);
}

void SubstitutionModel::organizeRateVectors(int seqLen, int numStates) {
  rateVectors.organize(seqLen, numStates);
}

RateVector* SubstitutionModel::selectRateVector(rv_request rq) {
  /*
   * This is a simple function right now but it will become hugely complex.
   * Given infomation about a BranchSegment and state of interest will return the corresponding rate vector.
   */
  return(rateVectors.select(rq));
}

// Getters

const double& SubstitutionModel::get_u() {
  return(u->getValue());
}

std::list<AbstractComponent*> SubstitutionModel::get_all_parameters() {
  std::list<AbstractComponent*> ret = {};
  for(auto it = rateVectors.col.begin(); it != rateVectors.col.end(); ++it) {
    for(auto jt = (*it)->rates.begin(); jt != (*it)->rates.end(); ++jt) {
      AbstractComponent* c = dynamic_cast<AbstractComponent*>(*jt);
      if(c == nullptr) {
	std::cerr << "Error: parameter not of AbstractComponent type." << std::endl;
	exit(EXIT_FAILURE);
      }
      ret.push_back(c);
    }
  }
  return(ret);
}

std::vector<RateVector*> SubstitutionModel::get_RateVectors() {
  return(rateVectors.col);
}

void SubstitutionModel::saveToFile(int gen, double l) {
  rateVectors.saveToFile(gen, l);
}

void SubstitutionModel::finalize() {
  // Add Indels to states table.
  states.state_to_int["-"] = -1;
  states.int_to_state[-1] = "-";

  rateVectors.Initialize();
}

// The ITERATOR

inline void SubstitutionModel::iterator::step_to_next_location() {
  ++location;
  while(location == location_iter_end) {
    endQ = step_to_next_component();
    if(endQ) {
      return;
    }
  }
}

inline bool SubstitutionModel::iterator::step_to_next_component() {
  cq.pop(); // Remove first element of queue.
  if(cq.empty()) {
    return(true);
  } else {
    // Add dependancies for new Abstract component at head of the queue.
    for(auto it = cq.front()->get_dependancies().begin(); it != cq.front()->get_dependancies().end(); ++it) {
      cq.push(*it);
    }
    Valuable* v = dynamic_cast<Valuable*>(cq.front());
    if(v != nullptr) {
      location = v->host_vectors.begin();
      location_iter_end = v->host_vectors.end();
    }
    return(false);
  }
}

SubstitutionModel::iterator::iterator(SubstitutionModel& s, bool e, AbstractComponent* modified_component) : sub_model(s), endQ(e) {
  cq = {};
  cq.push(modified_component);

  // This if statement is for cases when there are no parameters to be fitted.
  // Essentially a redundant if statement in vast majority of cases.
  if(cq.front() == nullptr) {
    endQ = true;
  } else {
    // Add dependancies to queue.
    for(auto it = cq.front()->get_dependancies().begin(); it != cq.front()->get_dependancies().end(); ++it) {
      cq.push(*it);
    }

    Valuable* v = dynamic_cast<Valuable*>(cq.front());
    if(v != nullptr) {
      location = v->host_vectors.begin();
      location_iter_end = v->host_vectors.end();
    }

    while(location == location_iter_end) {
      endQ = step_to_next_component();
      if(endQ) {
      return;
      }
    }
  }
}

SubstitutionModel::iterator& SubstitutionModel::iterator::operator++() {
  step_to_next_location();
  return(*this);
}

const rv_loc& SubstitutionModel::iterator::operator*() const {
  return(*location);
}

bool SubstitutionModel::iterator::at_end() const {
  return(endQ);
}

SubstitutionModel::iterator SubstitutionModel::modified_begin(AbstractComponent* modified_component) {
  return(SubstitutionModel::iterator(*this, false, modified_component));
}

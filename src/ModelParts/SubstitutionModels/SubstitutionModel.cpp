#include <iostream>

#include "SubstitutionModel.h"
#include "../../Environment.h"
#include "../../IO/Files.h"

extern Environment env;
extern IO::Files files;

SubstitutionModel::SubstitutionModel(Valuable* u) : u(u) {
  states.n = 0;
}

void SubstitutionModel::add_state(std::string s) {
  // Check if state is gap.
  states.possible.insert(s);
  states.state_to_int[s] = states.n;
  states.int_to_state[states.n] = s;
  states.n++;
}

RateVector* SubstitutionModel::create_rate_vector(IO::raw_rate_vector rv, Valuable* u) {
  std::vector<Valuable*> rates(states.n, nullptr);
  int s = states.state_to_int[rv.uc.state];

  // Setup Virtual Substitution rate.
  auto ptr_vir = rv.rates.begin();
  int i = 0;
  while(i != s) {
    ++ptr_vir;
    i++;
  }

  // Test code.
  VirtualSubstitutionRate* vir_rate = dynamic_cast<VirtualSubstitutionRate*>(*ptr_vir);

  if(vir_rate == nullptr) {
    std::cerr << "Error: expecting virtual substitution rate at position " << i << " in rate vector " << rv.name << "." << std::endl;
    exit(EXIT_FAILURE);
  }
  
  vir_rate->set_u(u);

  // Add each of the remaining parameters to the rate vector.
  for(int i = 0; i < states.n; i++) {
    AbstractComponent* param = rv.rates.front();
    if(param == nullptr) {
      std::cerr << "Error: nullptr for parameter in rate vector " << rv.name << " at position " << i << "." << std::endl;
      exit(EXIT_FAILURE);
    }

    if(i != s) {
      Valuable* v = dynamic_cast<Valuable*>(param);
      if(v == nullptr) {
	std::cerr << "Error: parameter in raw_rate_vector is not Valuable." << std::endl;
	exit(EXIT_FAILURE);
      }
      rates[i] = v;
      // Add dependancies.
      vir_rate->add_rate(rates[i]);
    } else {
      // Virtual Substitution rate.
      rates[i] = vir_rate;
    }
    rv.rates.pop_front();
  }

  if(not rv.rates.empty()) {
    std::cerr << "Error: unexpected number of rates - " << states.n << " is expected (given the number of states), however " << rv.rates.size() << " additional rate is found." << std::endl;
    exit(EXIT_FAILURE);
  }

  RateVector* new_rv = new RateVector(rv.name, rv.uc.state, &states, rates);
  return(new_rv);
}

void SubstitutionModel::from_raw_model(IO::raw_substitution_model* raw_sm) {
  // Read in states.
  for(auto it = raw_sm->states.begin(); it != raw_sm->states.end(); ++it) {
	add_state(*it);
  }

  // Read in rate vectors.
  for(auto raw_rv = raw_sm->rv_list.begin(); raw_rv != raw_sm->rv_list.end(); ++raw_rv) {
    RateVector* rv = create_rate_vector(*raw_rv, u);
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
  return(u->get_value());
}

void add_all_dependancies(std::list<AbstractComponent*>& all, AbstractComponent* parameter) {
  all.push_back(parameter);
  for(auto it = parameter->get_dependancies().begin(); it != parameter->get_dependancies().end(); ++it) {
    add_all_dependancies(all, *it);
  }
}

std::list<AbstractComponent*> SubstitutionModel::get_all_parameters() {
  std::list<AbstractComponent*> all_deps = {};
  for(auto it = rateVectors.col.begin(); it != rateVectors.col.end(); ++it) {
    for(auto jt = (*it)->rates.begin(); jt != (*it)->rates.end(); ++jt) {
      AbstractComponent* parameter = dynamic_cast<AbstractComponent*>(*jt);
      if(parameter == nullptr) {
	std::cerr << "Error: parameter not of AbstractComponent type." << std::endl;
	exit(EXIT_FAILURE);
      }
      add_all_dependancies(all_deps, parameter);
    }
  }
  return(all_deps);
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

  rateVectors.Initialize(&states);
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
  current_parameter++;
  if(current_parameter == valuables_end) { 
    //std::cout << std::endl;
    return(true);
  } else {
    location = (*current_parameter)->host_vectors.begin();
    location_iter_end = (*current_parameter)->host_vectors.end();
    return(false);
  }
}

SubstitutionModel::iterator::iterator(SubstitutionModel& s, bool e, AbstractComponent* modified_component) : sub_model(s), endQ(e) {
  current_parameter = modified_component->get_valuable_dependents().begin();
  valuables_end = modified_component->get_valuable_dependents().end();

  //AbstractComponent* comp = dynamic_cast<AbstractComponent*>(*current_parameter);
  //std::cout << "Begin: " << comp->get_name() << "-" << (*current_parameter)->host_vectors.size() << " ";

  //std::cout << "VD " << modified_component->get_valuable_dependents().size() << " : ";
  //for(auto it = modified_component->get_valuable_dependents().begin();
  //    it != modified_component->get_valuable_dependents().end(); ++it) {
  //  AbstractComponent* comp = dynamic_cast<AbstractComponent*>(*it);
  //  std::cout << comp->name << " ";
  //}
  //std::cout << std::endl;

  // This if statement is for cases when there are no parameters to be fitted.
  // Essentially a redundant if statement in vast majority of cases.
  if(current_parameter == valuables_end) {
    endQ = true;
  } else {
    location = (*current_parameter)->host_vectors.begin();
    location_iter_end = (*current_parameter)->host_vectors.end();

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
  //AbstractComponent* comp = dynamic_cast<AbstractComponent*>(*current_parameter);
  //std::cout << comp->name << " ";
  return(*location);
}

bool SubstitutionModel::iterator::at_end() const {
  return(endQ);
}

SubstitutionModel::iterator SubstitutionModel::modified_begin(AbstractComponent* modified_component) {
  return(SubstitutionModel::iterator(*this, false, modified_component));
}

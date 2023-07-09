#include <iostream>
#include <vector>

#include "SubstitutionModel.h"

#include "../../Environment.h"
#include "../../IO/Files.h"

#include "Parameters.h"

extern Environment env;
extern IO::Files files;

SubstitutionModel::SubstitutionModel(Valuable* u) : u(u) {
}

RateVector* SubstitutionModel::create_rate_vector(IO::raw_rate_vector rv, Valuable* u) {
  States* domain_states = &all_state_domains[rv.uc.domain];

  std::vector<Valuable*> rates(domain_states->n, nullptr);
  int s = domain_states->state_to_int[rv.uc.state];

  // Find and configure the Virtual Substitution rate.
  auto ptr_vir = rv.rates.begin();
  int i = 0;
  while(i != s) {
    ++ptr_vir;
    i++;
  }

  VirtualSubstitutionRate* vir_rate = dynamic_cast<VirtualSubstitutionRate*>(*ptr_vir);
  
  if(vir_rate == nullptr) {
    std::cerr << "Error: expecting virtual substitution rate at position " << i << " in rate vector " << rv.name << "." << std::endl;
    exit(EXIT_FAILURE);
  }
  
  vir_rate->set_u(u);

  // Add each of the remaining parameters to the rate vector.
  for(unsigned int i = 0; i < domain_states->n; i++) {
    AbstractComponent* param = rv.rates.front();
    if(param == nullptr) {
      std::cerr << "Error: nullptr for parameter in rate vector " << rv.name << " at position " << i << "." << std::endl;
      exit(EXIT_FAILURE);
    }

    if((int)i != s) {
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
    std::cerr << "Error: unexpected number of rates - " << domain_states->n << " is expected (given the number of states), however " << rv.rates.size() << " additional rate is found." << std::endl;
    exit(EXIT_FAILURE);
  }

  RateVector* new_rv = new RateVector(rv.name, rv.uc.state, domain_states, rates);
  return(new_rv);
}

void SubstitutionModel::configure_RateVectors(std::list<IO::raw_rate_vector> rv_list) {
  for(auto raw_rv = rv_list.begin(); raw_rv != rv_list.end(); ++raw_rv) {
    RateVector* rv = create_rate_vector(*raw_rv, u);
    rateVectors.add(rv, (*raw_rv).uc);
  }

  rateVectors.Initialize(all_state_domains);
}

void SubstitutionModel::configure_States(std::map<std::string, std::list<std::string>> raw_state_domains) {
  for(auto it = raw_state_domains.begin(); it != raw_state_domains.end(); ++it) {
    States new_state_domain = {};
    for(auto s = it->second.begin(); s != it->second.end(); ++s) {
      new_state_domain = add_to_States(new_state_domain, *s);
    }
    new_state_domain.state_to_int["-"] = -1;
    new_state_domain.int_to_state[-1] = "-";

    all_state_domains[it->first] = new_state_domain;
  }
}

void SubstitutionModel::from_raw_model(IO::raw_substitution_model* raw_sm) {
  configure_States(raw_sm->get_all_states());

  configure_RateVectors(raw_sm->get_rate_vector_list());

 // Parameter's counts out.
  files.add_file("parameters_counts_out", env.get<std::string>("OUTPUT.parameters_counts_out_file"), IOtype::OUTPUT);

  std::ostringstream counts_buffer;
  counts_buffer << "I,GEN,LogL";
  std::list<AbstractComponent*> all_parameters = get_all_parameters();
  for(auto it = all_parameters.begin(); it != all_parameters.end(); it++) {
    if((*it)->get_hidden() != true) {
      if(dynamic_cast<Valuable*>(*it) != nullptr) {
	counts_buffer << "," << (*it)->get_name();
      }
    }
  }
  counts_buffer << std::endl;

  files.write_to_file("parameters_counts_out", counts_buffer.str());

  //print_States(states);
}

const States* SubstitutionModel::get_state_domain(std::string domain_name) {
  auto s = all_state_domains.find(domain_name);
  if(s == all_state_domains.end()) {
    std::cerr << "Error: the domain \"" << domain_name << "\" is not recognized by the substitution model." << std::endl;
    exit(EXIT_FAILURE);
  }
  return(&(s->second));
}

std::map<std::string, States> SubstitutionModel::get_all_states() {
  return(all_state_domains);
}

void SubstitutionModel::organizeRateVectors() {
  rateVectors.organize();
}

RateVector* SubstitutionModel::selectRateVector(rv_request rq) {
  /*
   * This is a simple function right now but it will become hugely complex.
   * Given infomation about a BranchSegment and state of interest will return the corresponding rate vector.
   */
  return(rateVectors.select(rq));
}

unsigned long SubstitutionModel::get_hash_state(std::map<std::string, signed char> states) {
  return(rateVectors.get_hypothetical_hash_state(states));
}

// Getters

const double& SubstitutionModel::get_u() {
  return(u->get_value());
}

void add_all_dependancies(std::list<AbstractComponent*>& all, std::set<AbstractComponent*>& prev_parameters, AbstractComponent* parameter) {
  if(prev_parameters.find(parameter) == prev_parameters.end()) {
    prev_parameters.insert(parameter);
    all.push_back(parameter);
  }

  for(auto it = parameter->get_dependancies().begin(); it != parameter->get_dependancies().end(); ++it) {
    add_all_dependancies(all, prev_parameters, *it);
  }
}

std::list<AbstractComponent*> SubstitutionModel::get_all_parameters() {
  std::set<AbstractComponent*> prev_parameters = {};
  std::list<AbstractComponent*> all_deps = {};
  for(auto it = rateVectors.col.begin(); it != rateVectors.col.end(); ++it) {
    for(auto jt = (*it)->rates.begin(); jt != (*it)->rates.end(); ++jt) {
      AbstractComponent* parameter = dynamic_cast<AbstractComponent*>(*jt);
      if(parameter == nullptr) {
	std::cerr << "Error: parameter not of AbstractComponent type." << std::endl;
	exit(EXIT_FAILURE);
      }
      add_all_dependancies(all_deps, prev_parameters, parameter);
    }
  }
  return(all_deps);
}

std::vector<RateVector*> SubstitutionModel::get_RateVectors() {
  return(rateVectors.col);
}

void SubstitutionModel::saveToFile(uint128_t gen, double l, std::map<RateVector*, std::vector<double>> counts_by_rv) {
  rateVectors.saveToFile(gen, l);

  static int i = -1;
  ++i;
  
  // Parameter's substitution counts.
  std::string line = std::to_string(i) + "," + gen.str() + "," + std::to_string(l);
  std::list<AbstractComponent*> all_parameters = get_all_parameters();
  for(auto it = all_parameters.begin(); it != all_parameters.end(); it++) {
    if((*it)->get_hidden() != true) {
      Valuable* v = dynamic_cast<Valuable*>(*it);
      if(v != nullptr) {
	double total = 0.0;
	const std::list<rv_loc> host_rvs = rateVectors.get_host_vectors(v);
	for(auto it = host_rvs.begin(); it != host_rvs.end(); ++it) {
	  total += counts_by_rv[it->rv][it->pos];
	}
	line += "," + std::to_string(total);
      }
    }
  }

  files.write_to_file("parameters_counts_out", line + "\n");
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
    location = sub_model.rateVectors.get_host_vectors(*current_parameter).begin();
    location_iter_end = sub_model.rateVectors.get_host_vectors(*current_parameter).end();
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
    location = sub_model.rateVectors.get_host_vectors(*current_parameter).begin();
    location_iter_end = sub_model.rateVectors.get_host_vectors(*current_parameter).end();
    
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

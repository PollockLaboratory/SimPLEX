#include "RateVector.h"

#include "../../Environment.h"
#include "../../IO/Files.h"

extern Environment env;
extern IO::Files files;

int RateVector::IDc = 0;

// Figure out locations Data Structure.
RateVector::RateVector(std::string name, std::string s, const States* states, std::vector<Valuable*> params) : name(name), states(states) {
  state = states->state_to_int.at(s);
  id = IDc++;
  rates = params;
}

float RateVector::operator[](int i) {
  
  return(rates[i]->get_value());
}

float RateVector::get_rate_ratio(int i) {
  return(rates[i]->get_value() / rates[i]->get_old_value());
}

const int& RateVector::getID() {
  return(id);
}

int RateVector::size() {
  return(rates.size());
}

const std::string& RateVector::get_name() {
  return(name);
}

std::string RateVector::get_state() {
  return(states->int_to_state.at(state));
}

std::string RateVector::get_state_by_pos(int pos) {
  return(states->int_to_state.at(pos));
}

// Util.
void RateVector::print() {
  std::cout << "RateVector:\t" << name << "\t";
  for(auto it = rates.begin(); it != rates.end(); ++it) {
    std::cout << (*it)->get_value() << " ";
  }
  std::cout << std::endl;
}

ExtendedState ExtendedState_Null() {
  return(std::map<std::string, std::string>({}));
}

std::string ExtendedState_toString(ExtendedState state) {
  std::string ret = "";
  for(auto it = state.begin(); it != state.end(); ++it) {
    if(ret == "") {
      ret = it->second;
    } else {
      ret += "+" + it->second;
    }
  }
  return(ret);
}

// COLLECTIONS of rate vectors.
RateVectorSet::RateVectorSet() {
  n_domains = 0;
  domain_names = {};
}

void RateVectorSet::Initialize(States states, std::map<std::string, States> hidden_states) {
  // Set States.
  all_states = hidden_states;

  // Set Domain names.
  for(auto it = all_states.begin(); it != all_states.end(); ++it) {
    domain_names.push_back(it->first);
  }
  
  // Output file.
  files.add_file("rate_vectors_out", env.get<std::string>("OUTPUT.rate_vectors_out_file"), IOtype::OUTPUT);

  std::ostringstream buffer;
  buffer << "I,GEN,LogL,NAME,ANC";
  // the count is reduced by -1 to take into account gaps.
  for(unsigned int i = 0; i < states.int_to_state.size() - 1; ++i) {
    buffer << "," << states.int_to_state[i];
  }
  buffer << std::endl;

  files.write_to_file("rate_vectors_out", buffer.str());
}

void RateVectorSet::add(RateVector* rv, IO::rv_use_class uc) {
  // Set the pointer inside the AbstractValue to point to the host rate vector.
  for(unsigned int i = 0; i < rv->rates.size(); i++) {
    // Set up the AbstractValues themselves.
    //rv->rates[i]->add_host_vector(rv, i);
    Valuable* v = rv->rates[i];
    parameter_locations[v].push_back({rv, i});
  }

  // Add rate vector to base collection.
  col.push_back(rv);

  // Add RateVector to collection in RateVectorSet.
  id_to_uc[rv->getID()] = uc;
}

RateVector* RateVectorSet::select(rv_request rq) {
  if(rq.states.empty() == true) {
    std::cerr << "Sort this out." << std::endl;
    exit(EXIT_FAILURE);
  }

  ExtendedState state = ExtendedState_Null();
  for(auto it = rq.states.begin(); it != rq.states.end(); ++it) {
    state[it->first] = all_states[it->first].int_to_state[it->second];
  }

  //std::cout << "Select: " << rq.state << " " << rq.domain << " " << ExtendedState_toString(state) << std::endl;
  RateVector* rv = ex_state_to_rv[rq.domain][state];
  if(rv == nullptr) {
    std::cout << "Error: attempting to dispatch nullptr instead of RateVector*" << std::endl;
    exit(EXIT_FAILURE);
  }
  return(rv);
}

std::list<ExtendedState> expandStates(std::list<ExtendedState> base, std::map<std::string, States> extension) {
  std::list<ExtendedState> new_states = base;
  for(auto domain_it = extension.begin(); domain_it != extension.end(); ++domain_it) {
    std::list<ExtendedState> extended_states = {};
    for(auto it = new_states.begin(); it != new_states.end(); ++it) {
      for(auto jt = domain_it->second.possible.begin(); jt != domain_it->second.possible.end(); ++jt) {
	ExtendedState s = *it;
	s[domain_it->first] = *jt;
	extended_states.push_back(s);
      }
    }
    new_states = extended_states;
  }
  return(new_states);
}

void RateVectorSet::organize(int seqLen) {
  /*
   * Organizes RateVectors into logical structure now tree data is known.
   */

  // Looks at possible states.
  //std::cout << "Domains: [ ";
  //for(auto it = domain_names.begin(); it != domain_names.end(); ++it) {
  // std::cout << *it << " ";
  //}
  //std::cout << "]" << std::endl;

  // Calculated extended states.
  std::list<ExtendedState> extended_states = expandStates({ExtendedState_Null()}, all_states);

  //std::cout << "Extended States: [ ";
  //for(auto it = extended_states.begin(); it != extended_states.end(); ++it) {
  // std::cout << ExtendedState_toString(*it) << " ";
  //}
  //std::cout << "]" << std::endl;

  // Make empty tree.
  for(auto domain = domain_names.begin(); domain != domain_names.end(); ++domain) {
    ex_state_to_rv[*domain] = {};
    for(auto state = extended_states.begin(); state != extended_states.end(); ++state) {
      ex_state_to_rv[*domain][*state] = nullptr;
    }
  }

  // Fill tree - loop through vector.
  for(auto it = col.begin(); it != col.end(); ++it) {
    if(*it == nullptr) {
      std::cout << "Error: RateVector is nullptr." << std::endl;
      exit(EXIT_FAILURE);
    }
    IO::rv_use_class uc = id_to_uc[(*it)->getID()];
    ExtendedState s = ExtendedState_Null();
    s[uc.domain] = uc.state;
    std::list<ExtendedState> applicable_states = {s};

    // Find all Extended states that apply.
    for(auto it = uc.secondary_state.begin(); it != uc.secondary_state.end(); ++it) {
      if(it->second == "*") {
	// Applies to all states in domain.
	std::map<std::string, States> extension_domain = {};
	extension_domain[it->first] = all_states[it->first];
	applicable_states = expandStates(applicable_states, extension_domain);
      } else {
	// Applies to specific states in domain.
	for(auto jt = applicable_states.begin(); jt != applicable_states.end(); ++jt) {
	  (*jt)[it->first] = it->second;
	}
      }
    }

    // Print for Debug.
    //std::cout << (*it)->get_name() << ": [ ";
    //for(auto jt = applicable_states.begin(); jt != applicable_states.end(); ++jt) {
    // std::cout << ExtendedState_toString(*jt) << " ";
    //}
    //std::cout << "]" << std::endl;

    // Set ptrs to rate vectors based on applicable state.
    for(auto jt = applicable_states.begin(); jt != applicable_states.end(); ++jt) {
      if(ex_state_to_rv[uc.domain][*jt] != nullptr) {
	std::cerr << "Error: rate vector has already been assigned." << std::endl;
	exit(EXIT_FAILURE);
      } else {
	ex_state_to_rv[uc.domain][*jt] = *it;
      }
    }
  }

  // Check there are no nullpts remaining.
  for(auto domain = domain_names.begin(); domain != domain_names.end(); ++domain) {
    //std::cout << "Domain check: " << *domain << " " << extended_states.empty() << std::endl;
    for(auto state = extended_states.begin(); state != extended_states.end(); ++state) {
      //std::cout << "Check: " << *domain << " " << ExtendedState_toString(*state) << " " << (ex_state_to_rv[*domain][*state] == nullptr) << std::endl;
      if(ex_state_to_rv[*domain][*state] == nullptr) {
	std::cerr << "Error: no rate vector specified for " << ExtendedState_toString(*state) << " in domain \"" << *domain << "\"." << std::endl;
	exit(EXIT_FAILURE);
      }
    }
  } 
}

const std::list<rv_loc>& RateVectorSet::get_host_vectors(Valuable* v) {
  return(parameter_locations[v]);
}

void RateVectorSet::print() {
  for(std::vector<RateVector*>::iterator it = col.begin(); it != col.end(); ++it) {
    (*it)->print();
  }
}

void RateVectorSet::saveToFile(int gen, double l) {
  static int i = -1;
  ++i;

  std::ostringstream buffer;
  for(auto it = col.begin(); it != col.end(); ++it) {
    buffer << i << "," << gen << "," << l << "," << (*it)->get_name() << "," << env.integer_to_state[(*it)->state];
    for(auto jt = (*it)->rates.begin(); jt != (*it)->rates.end(); ++jt) {
      buffer << "," << (*jt)->get_value();
    }
    buffer << std::endl;
  }
  files.write_to_file("rate_vectors_out", buffer.str());
}

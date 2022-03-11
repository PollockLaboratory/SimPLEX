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

void RateVectorSet::Initialize(std::map<std::string, States> all_states) {
  // Set States.
  this->all_states = all_states;

  // Set Domain names.
  for(auto it = all_states.begin(); it != all_states.end(); ++it) {
    domain_names.push_back(it->first);
  }
  
  // Output file.
  files.add_file("rate_vectors_out", env.get<std::string>("OUTPUT.rate_vectors_out_file"), IOtype::OUTPUT);

  std::ostringstream buffer;
  buffer << "I,GEN,LogL,NAME,ANC,VECTOR" << std::endl;
  // the count is reduced by -1 to take into account gaps.

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
  RateVector* rv = state_to_rv[rq.domain][rq.ex_state];

  if(rv == nullptr) {
    std::cerr << "Error: attempting to dispatch nullptr instead of RateVector*" << std::endl;
    exit(EXIT_FAILURE);
  }

  return(rv);
}

unsigned long hash(std::list<signed char> seq) {
  unsigned long hash = 5381;
  int c;

  for(auto it = seq.begin(); it != seq.end(); ++it) {
    c = *it;
    hash = ((hash << 5) + hash) + c; /* hash * 33 + c */
  }

  return(hash);
}

unsigned long hash_step(unsigned long h, signed char c) {
  return(((h << 5) + h) + c);
}

unsigned long RateVectorSet::get_hash_state(const std::map<std::string, std::vector<signed char>*>& sequences, int pos) const {
  unsigned long h = 5381;
  for(auto domain_it = all_states.begin(); domain_it != all_states.end(); ++domain_it) {
    const std::vector<signed char>* s = sequences.at(domain_it->first);
    h = hash_step(h, (*s)[pos]);
  }

  return(h);
}

unsigned long RateVectorSet::get_hypothetical_hash_state(const std::map<std::string, std::vector<signed char>*>& sequences, int pos, std::string domain_name, signed char state) const {
  unsigned long h = 5381;
  signed char c;
  for(auto domain_it = all_states.begin(); domain_it != all_states.end(); ++domain_it) {
    if(domain_it->first == domain_name) {
      c = state;
    } else {
      const std::vector<signed char>* s = sequences.at(domain_it->first);
      c = (*s)[pos];
    }
    h = hash_step(h, c);
  }

  return(h);
}

std::list<signed char> RateVectorSet::ex_to_list(ExtendedState ex) {
  std::list<signed char> l = {};
  for(auto domain_it = all_states.begin(); domain_it != all_states.end(); ++domain_it) {
    l.push_back(domain_it->second.state_to_int[ex[domain_it->first]]);
  }
  return(l);
}

std::list<std::list<signed char>> RateVectorSet::configure_hash(std::map<std::string, States> all_states) {
  std::set<unsigned long> hash_values = {};

  std::list<std::list<signed char>> new_states = {{}};
  for(auto domain_it = all_states.begin(); domain_it != all_states.end(); ++domain_it) {
    std::list<std::list<signed char>> tmp_states = {};
    for(auto it = new_states.begin(); it != new_states.end(); ++it) {
      for(auto jt = domain_it->second.possible.begin(); jt != domain_it->second.possible.end(); ++jt) {
	std::list<signed char> s = *it;
	s.push_back(domain_it->second.state_to_int[*jt]);

	// Check for hash collision.
	if(hash_values.find(hash(s)) != hash_values.end()) {
	  std::cerr << "Error: collision in hash values.\n" << std::endl;
	  exit(EXIT_FAILURE);
	}

	tmp_states.push_back(s);
      }
    }
    new_states = tmp_states;
  }

  return(new_states);
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

void RateVectorSet::organize() {
  /*
   * Organizes RateVectors into logical structure now tree data is known.
   */

  std::list<std::list<signed char>> tmp_states = configure_hash(all_states);

  // Make empty tree
  for(auto domain = domain_names.begin(); domain != domain_names.end(); ++domain) {
    state_to_rv[*domain] = {};
    for(auto state = tmp_states.begin(); state != tmp_states.end(); ++state) {
      state_to_rv[*domain][hash(*state)] = nullptr;
    }
  }

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

    // Set ptrs to rate vectors based on applicable state.
    for(auto jt = applicable_states.begin(); jt != applicable_states.end(); ++jt) {
      unsigned long h = hash(ex_to_list(*jt));
      if(state_to_rv[uc.domain][h] != nullptr) {
	std::cerr << "Error: rate vector has already been assigned." << std::endl;
	exit(EXIT_FAILURE);
      } else {
	state_to_rv[uc.domain][h] = *it;
      }
    }
  }

  // Check there are no nullpts remaining.
  for(auto domain = domain_names.begin(); domain != domain_names.end(); ++domain) {
    //std::cout << "Domain check: " << *domain << " " << extended_states.empty() << std::endl;
    for(auto state = tmp_states.begin(); state != tmp_states.end(); ++state) {
      //std::cout << "Check: " << *domain << " " << ExtendedState_toString(*state) << " " << (ex_state_to_rv[*domain][*state] == nullptr) << std::endl;
      if(state_to_rv[*domain][hash(*state)] == nullptr) {
	std::cerr << "Error: no rate vector specified for " << hash(*state) << " in domain \"" << *domain << "\"." << std::endl;
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

void RateVectorSet::saveToFile(uint128_t gen, double l) {
  static int i = -1;
  ++i;

  std::ostringstream buffer;
  for(auto it = col.begin(); it != col.end(); ++it) {
    buffer << i << "," << gen << "," << l << "," << (*it)->get_name() << "," << (*it)->get_state();
    for(auto jt = (*it)->rates.begin(); jt != (*it)->rates.end(); ++jt) {
      buffer << "," << (*jt)->get_value();
    }
    buffer << std::endl;
  }
  files.write_to_file("rate_vectors_out", buffer.str());
}

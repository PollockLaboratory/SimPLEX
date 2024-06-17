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
  static_domains = {};
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

void RateVectorSet::add(RateVector* rv, IO::RVScope scope) {
  // Set the pointer inside the AbstractValue to point to the host rate vector.
  for(unsigned int i = 0; i < rv->rates.size(); i++) {
    // Set up the AbstractValues themselves.
    //rv->rates[i]->add_host_vector(rv, i);
    Valuable* v = rv->rates[i];
    parameter_locations[v].push_back({rv, (int)i});
  }

  // Add rate vector to base collection.
  collection.push_back(rv);

  // Add RateVector to collection in RateVectorSet.
  id_to_uc[rv->getID()] = scope;
}

void RateVectorSet::mark_static_state(std::string state_domain) {
  // This marks a state as static/unsampled and therefore there is no need for rate vectors.
  this->static_domains.insert(state_domain);
}

bool RateVectorSet::is_static(std::string state_domain) {
  return(this->static_domains.find(state_domain) != this->static_domains.end());
}

RateVector* RateVectorSet::select(RVQuery query) {
  RateVector* rv = state_to_rv[query.domain][query.ex_state];

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
    //std::cout << domain_it->first << " ";
    if(domain_it->first == domain_name) {
      c = state;
      //std::cout << "primary " << (unsigned int)c << std::endl;
    } else {
      const std::vector<signed char>* s = sequences.at(domain_it->first);
      c = (*s)[pos];
      //std::cout << "alt " << (unsigned int)c << std::endl;
    }
    h = hash_step(h, c);
  }

  return(h);
}


unsigned long RateVectorSet::get_hypothetical_hash_state(std::map<std::string, signed char> states) {
  unsigned long h = 5381;
  signed char c;

  for(auto domain_it = all_states.begin(); domain_it != all_states.end(); ++domain_it) {
    c = states[domain_it->first];
    h = hash_step(h, c);
  }
  return(h);
}

std::list<signed char> RateVectorSet::flatten_compound_state(ExtendedState compound_state) {
  /*
   * Uses the order of this->all_states to flatten ExtendedState.
   */
  std::list<signed char> l = {};
  for(auto& [state_domain, states] : this->all_states) {
    l.push_back(states.state_to_int[compound_state[state_domain]]);
  }
  return(l);
}

std::list<std::list<signed char>> enumerate_state_combinations(std::map<std::string, States> all_states) {
  /*
   * Ennumerates all the possible state combinations.
   * e.g. two states NUCLEOTIDE [ A T C G ] and SECONDARY [ X Y ]
   * the list { { A, X } { A, Y }, { T, X, }, { T, Y } ... } will be generated.
   * however rather than state name being in the list such as A, the corresponding index will be used.
   * e.g. { { 0, 0 }, { 0, 1 }, { 1, 0 }, { 1, 1 }, etc
   */
  
  std::list<std::list<signed char>> acc_states = {{}};
  std::set<unsigned long> hash_values = {};

  for (auto& [domain_name, states] : all_states) {
    std::list<std::list<signed char>> intermediate_acc = {};

    for(std::list<signed char> acc : acc_states) {
      for (const std::string& state : states.possible) {
        std::list<signed char> new_state = acc;
        new_state.push_back(states.state_to_int[state]);

        // Check for hash collisions - this shouldn't be able to happen but to be safe.
        if(hash_values.find(hash(new_state)) != hash_values.end()) {
          std::cerr << "Error: collision in hash values.\n" << std::endl;
          exit(EXIT_FAILURE);
        }

        hash_values.insert(hash(new_state));
        intermediate_acc.push_back(new_state);
      }
    }

    acc_states = intermediate_acc;
  }

  return(acc_states);
}

std::list<ExtendedState> expandStates(std::list<ExtendedState> base, std::map<std::string, States> extension) {
  /*
   * expands a compound state across a full state domain. e.g.
   * takes { { A }, { B } } and extends across new state domain { X, Y }
   * results in { { A, X }, { A, Y }, { B, X }, { B, Y } }
   */
  std::list<ExtendedState> possible_compound_states = base;
  for(const auto& [ state_domain, states ] : extension) {
    std::list<ExtendedState> extended_states = {};
    for(const auto& acc_state : possible_compound_states) {
      for(auto jt : states.possible) {
        ExtendedState new_compound_state = acc_state;
        new_compound_state[state_domain] = jt;
        extended_states.push_back(new_compound_state);
      }
    }
    possible_compound_states = extended_states;
  }

  return(possible_compound_states);
}

void RateVectorSet::organize() {
  /*
   * Organizes RateVectors into logical structure now tree data is known.
   */

  std::list<std::list<signed char>> all_possible_states = enumerate_state_combinations(all_states);

  // Make empty dictionary/map
  // state_domain->hash(state)->RateVector*
  for (std::string& domain : domain_names) {
    this->state_to_rv[domain] = {};
    for (std::list<signed char>& state : all_possible_states) this->state_to_rv[domain][hash(state)] = nullptr;
  }

  for (RateVector* rate_vector : this->collection) {
    if(rate_vector == nullptr) {
      std::cout << "Error: RateVector is nullptr." << std::endl;
      exit(EXIT_FAILURE);
    }

    IO::RVScope scope = this->id_to_uc[rate_vector->getID()];
    std::map<std::string, std::string> s({});
    s[scope.domain] = scope.state;
    std::list<ExtendedState> applicable_states = { s };

    // Find all Extended states that apply.
    for(const auto& [state_domain, valid_context] : scope.secondary_state) {
      if(valid_context == "*") {
        // Applies to all states in domain.
        std::map<std::string, States> extension_domain = {};
        extension_domain[state_domain] = all_states[state_domain];
        applicable_states = expandStates(applicable_states, extension_domain);
      } else {
        // Applies to specific states in domain.
        for(auto& jt : applicable_states) {
          jt[state_domain] = valid_context;
        }
      }
    }

    // Set ptrs to rate vectors based on applicable state.
    for(auto& compund_state : applicable_states) {
      unsigned long h = hash(flatten_compound_state(compund_state));
      if(state_to_rv[scope.domain][h] != nullptr) {
        std::cerr << "Error: rate vector has already been assigned." << std::endl;
        exit(EXIT_FAILURE);
      } else {
        state_to_rv[scope.domain][h] = rate_vector;
      }
    }
  }

  // Check there are no nullpts remaining.
  for(const auto& domain : domain_names) {
    if (this->is_static(domain)) continue;
    for(const auto& state : all_possible_states) {
      if(state_to_rv[domain][hash(state)] == nullptr) {
        std::cerr << "Error: no rate vector specified for " << hash(state) << " in domain \"" << domain << "\"." << std::endl;
        exit(EXIT_FAILURE);
      }
    }
  } 
}

const std::list<rv_loc>& RateVectorSet::get_host_vectors(Valuable* v) {
  return(parameter_locations[v]);
}

void RateVectorSet::print() {
  for(std::vector<RateVector*>::iterator it = collection.begin(); it != collection.end(); ++it) {
    (*it)->print();
  }
}

void RateVectorSet::saveToFile(uint128_t gen, double l) {
  static int i = -1;
  ++i;

  std::ostringstream buffer;
  for(auto it = collection.begin(); it != collection.end(); ++it) {
    buffer << i << "," << gen << "," << l << "," << (*it)->get_name() << "," << (*it)->get_state();
    for(auto jt = (*it)->rates.begin(); jt != (*it)->rates.end(); ++jt) {
      buffer << "," << (*jt)->get_value();
    }
    buffer << std::endl;
  }
  files.write_to_file("rate_vectors_out", buffer.str());
}

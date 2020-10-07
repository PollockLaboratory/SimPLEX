#include "RateVector.h"
#include "../Sequence.h"
#include "../Trees/TreeParts.h"
#include "../../Environment.h"
#include "../../IO/Files.h"

#include <sstream>

extern Environment env;
extern IO::Files files;

int RateVector::IDc = 0;

// Figure out locations Data Structure.
RateVector::RateVector(std::string name, std::string s, const States* states, std::vector<Valuable*> params) : name(name), states(states) {
  state = states->state_to_int.at(s);
  id = IDc++;
  //std::cout << id << " " << name << std::endl;
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

// COLLECTIONS of rate vectors.
RateVectorSet::RateVectorSet() {
}

void RateVectorSet::Initialize(States* states) {
  // Output file.
  files.add_file("rate_vectors_out", env.get<std::string>("OUTPUT.rate_vectors_out_file"), IOtype::OUTPUT);

  std::ostringstream buffer;
  buffer << "I,GEN,LogL,NAME,ANC";
  // the count is reduced by -1 to take into account gaps.
  for(int i = 0; i < states->int_to_state.size() - 1; ++i) {
    buffer << "," << states->int_to_state[i];
  }
  buffer << std::endl;

  files.write_to_file("rate_vectors_out", buffer.str());
}

void RateVectorSet::add(RateVector* rv, IO::rv_use_class uc) {
  // Set the pointer inside the AbstractValue to point to the host rate vector.
  for(unsigned int i = 0; i < rv->rates.size(); i++) {
    // Set up the AbstractValues themselves.
    rv->rates[i]->add_host_vector(rv, i);
  }
  // Add RateVector to collection in RateVectorSet.
  id_to_uc[rv->getID()] = uc;
  col.push_back(rv);
}

RateVector* RateVectorSet::select(rv_request rq) {
  std::list<RateVector*> possible_rv = rv_tree[rq.pos][rq.state];
  if(possible_rv.size() != 1) {
    std::cerr << "Error: Incorrect number of possible rate vectors." << std::endl;
    exit(EXIT_FAILURE);
  } else {
    RateVector* rv = possible_rv.front();
    if(rv->state != rq.state) {
      std::cerr << "Error: Trying to dispatch incorrect rate vector." << std::endl;
      exit(EXIT_FAILURE); 
    }
    return(possible_rv.front());
  }
}

void RateVectorSet::organize(int seqLen, int numStates) {
  /*
   * Organizes RateVectors into logical structure now tree data is known.
   */

  // Establish structure.
  rv_tree = std::vector<std::vector<std::list<RateVector*>>>(seqLen, std::vector<std::list<RateVector*>>(numStates));
  for(int l = 0; l < seqLen; l++) {
    for(int s = 0; s < numStates; s++) {
      rv_tree[l][s] = {};
    }
  }

  // Fill with pointers to RateVectors.
  for(auto it = col.begin(); it != col.end(); ++it) {
    // Add rate vector.
    IO::rv_use_class uc = id_to_uc[(*it)->getID()];
    // If no positions given then will apply to all positions.
    if(uc.pos.size() == 0) {
      for(int i = 0; i < seqLen; i++) {
	rv_tree[i][(*it)->state].push_back(*it);
      }
    } else {
      for(auto p = uc.pos.begin(); p != uc.pos.end(); ++p) {
	rv_tree[*p][(*it)->state].push_back(*it);
      }
    }
  }
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

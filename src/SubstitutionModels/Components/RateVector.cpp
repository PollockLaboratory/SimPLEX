#include "RateVector.h"
#include "../../Sequence.h"
#include "../../Trees/TreeParts.h"
#include "../../Environment.h"
#include "../../IO.h"

extern Environment env;
extern IO::Files files;

int RateVector::IDc = 0;

// Figure out locations Data Structure.
RateVector::RateVector(std::string name, int state, std::vector<AbstractValue*> params) : name(name), state(state) {
  size = params.size();
  rates = params;	
  //locations = {};
  counts = std::vector<int>(env.num_states, 0);
  logLikelihoods = std::vector<double>(env.num_states, 0);

  for(unsigned int i = 0; i < rates.size(); i++) {
    valueID_to_state[rates[i]->get_ID()] = i;
  }
}

float RateVector::operator[](int i) {
  
  return(rates[i]->getValue());
}

float RateVector::get_rate_ratio(int i) {
  return(rates[i]->getValue() / rates[i]->getOldValue());
}

// Util.
void RateVector::print() {
  std::cout << "RateVector:\t" << name << "\t";
  for(auto it = rates.begin(); it != rates.end(); ++it) {
    std::cout << (*it)->getValue() << " ";
  }
  std::cout << std::endl;
}

// COLLECTIONS of rate vectors.
std::ofstream RateVectorSet::out_file;

RateVectorSet::RateVectorSet() {
}

void RateVectorSet::Initialize() {
  files.add_file("rate_vectors", env.get<std::string>("OUTPUT.rate_vectors_out_file"), IOtype::OUTPUT);
  out_file = files.get_ofstream("rate_vectors");
		
  out_file << "I,GEN,LogL,NAME,ANC";
  for(auto it = env.state_to_integer.begin(); it != env.state_to_integer.end(); ++it) {
    out_file << "," << it->first;
  }
  out_file << std::endl;
}

RateVector*& RateVectorSet::operator[] (const int i) {
  if(col[i]->state != i) {
    std::cerr << "Error: RateVectorSet dispatching incorrect rate vector. " << std::endl;
    exit(EXIT_FAILURE);
  }
  return(col[i]);
}

void RateVectorSet::add(RateVector* rv) {
  // Set the pointer inside the AbstractValue to point to the host rate vector.
  for(unsigned int i = 0; i < rv->rates.size(); i++) {
    // Set up the AbstractValues themselves.
    rv->rates[i]->add_host_vector(rv, i);
  }
  // Add RateVector to collection in RateVectorSet. 
  col.push_back(rv);
}

void RateVectorSet::print() {
  for(std::vector<RateVector*>::iterator it = col.begin(); it != col.end(); ++it) {
    (*it)->print();
  }
}

void RateVectorSet::saveToFile(int gen, double l) {
  static int i = -1;
  ++i;
  for(auto it = col.begin(); it != col.end(); ++it) {
    out_file << i << "," << gen << "," << l << "," << (*it)->name << "," << (*it)->state;
    for(auto jt = (*it)->rates.begin(); jt != (*it)->rates.end(); ++jt) {
      out_file << "," << (*jt)->getValue();
    }
    out_file << std::endl;
  }
}

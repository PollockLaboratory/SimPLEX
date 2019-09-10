#include "RateVector.h"
#include "../../Sequence.h"
#include "../../Trees/TreeParts.h"
#include "../../Environment.h"
#include "../../IO/Files.h"

extern Environment env;
extern IO::Files files;

int RateVector::IDc = 0;

// Figure out locations Data Structure.
RateVector::RateVector(std::string name, int state, std::vector<AbstractValue*> params) : name(name), state(state) {
  id = IDc;
  IDc++;
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

const int& RateVector::getID() {
  return(id);
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
  std::cout << "Initializing rate vector" << std::endl;

  // Output file.
  files.add_file("rate_vectors", env.get<std::string>("OUTPUT.rate_vectors_out_file"), IOtype::OUTPUT);

  std::cout << "1." << std::endl;
  out_file = files.get_ofstream("rate_vectors");
		
  std::cout << "2." << std::endl;
  out_file << "I,GEN,LogL,NAME,ANC";
  for(auto it = env.state_to_integer.begin(); it != env.state_to_integer.end(); ++it) {
    out_file << "," << it->first;
  }
  out_file << std::endl;

  // Collection of Rate Vectors.

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
  std::cout << seqLen << " " << numStates << std::endl;
  for(int l = 0; l < seqLen; l++) {
    for(int s = 0; s < numStates; s++) {
      //std::cout << l << " " << s << std::endl;
      rv_tree[l][s] = {};
    }
  }

  // Fill with pointers to RateVectors.
  for(auto it = col.begin(); it != col.end(); ++it) {
    // Add rate vector.
    IO::rv_use_class uc = id_to_uc[(*it)->getID()];
    std::cout << "RateVector: " << (*it)->getID() << " " << uc.pos.size() << std::endl;
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
  for(auto it = col.begin(); it != col.end(); ++it) {
    out_file << i << "," << gen << "," << l << "," << (*it)->name << "," << (*it)->state;
    for(auto jt = (*it)->rates.begin(); jt != (*it)->rates.end(); ++jt) {
      out_file << "," << (*jt)->getValue();
    }
    out_file << std::endl;
  }
}

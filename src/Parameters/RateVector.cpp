#include "RateVector.h"
#include "Sequence.h"
#include "TreeParts.h"
#include "Environment.h"
#include "IO.h"

extern Environment env;
extern IO::Files files;

int RateVector::IDc = 0;

// Figure out locations Data Structure.
RateVector::RateVector(std::string name, int state, std::vector<AbstractValue*> params) : name(name), state(state) {
  size = params.size();
  rates = params;	
  locations = {};
  // TODO 20 should not be hard coded.
  counts = std::vector<int>(20, 0);
}

void RateVector::add_location(int pos, BranchSegment* bs) {
  bpos loc = {bs, pos};
  locations.insert(loc);
}

void RateVector::remove_location(int pos, BranchSegment* bs) {
  bpos loc = {bs, pos};
  //std::cout << "B: " << bs << " P: " << pos << std::endl;
  //for(auto l = locations.begin(); l != locations.end(); ++l) {
  // std::cout << "[ B: " << bs << " P: " << pos << " ] " << std::endl;
  // }
  int s = locations.size();
  locations.erase(loc);
  if(s == locations.size()) {
    std::cerr << "Error: location not removed from RateVector.locations." << std::endl;
    exit(EXIT_FAILURE);
  }
}

void RateVector::clear_locations() {
  locations = {};
}

std::unordered_set<bpos> RateVector::get_locations() {
  return(locations);
}

void RateVector::update_counts() {
  // Gets the counts of substatutions that this rate vector applies to.
  bpos b;
  int dec_state;
  substitution sub;
  // TODO 20 should not be hard coded.
  counts = std::vector<int>(20, 0);
  for(auto it = locations.begin(); it != locations.end(); ++it) {
    b = *it;
    dec_state = (*(b.branch->decendant->sequence))[b.pos];
    if(state != dec_state) {
      // Normal substitution - careful of indels.
      counts[dec_state]++;
    } else {
      // Check for virtual substitution.
      sub = b.branch->subs[b.pos];
      if(sub.pos != -1) {
	// -1 in pos indicates no substitution.
	counts[state]++;
      }
    }
  }
}

double RateVector::get_logLikelihood() {
  double logL = 0.0;
  for(int i = 0; i < rates.size(); i++) {
    logL += counts[i] * log(rates[i]->getValue());
  }
  return(logL);
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
  files.add_file("rate_vectors", env.get("rate_vectors_out_file"), IOtype::OUTPUT);
  out_file = files.get_ofstream("rate_vectors");
		
  out_file << "I,GEN,LogL,NAME,ANC";
  for(auto it = env.state_to_integer.begin(); it != env.state_to_integer.end(); ++it) {
    out_file << "," << it->first;
  }
  out_file << std::endl;
}

RateVector*& RateVectorSet::operator[] (const int i) {
  RateVector* rv = c[i];
  if(c[i]->state != i) {
    std::cerr << "Error: RateVectorSet dispatching incorrect rate vector. " << std::endl;
    exit(EXIT_FAILURE);
  }
  return(c[i]);
}

void RateVectorSet::add(RateVector* v) {
  for(int i = 0; i < v->rates.size(); i++) {
    v->rates[i]->state = i;
    v->rates[i]->rv = v;
  }
  c.push_back(v);
}

void RateVectorSet::get_counts() {
  // Calls get_counts in each rate vector.
  for(auto r = c.begin(); r != c.end(); ++r) {
    (*r)->update_counts();
  }
}

void RateVectorSet::clear_locations() {
  for(auto r = c.begin(); r != c.end(); ++r) {
    (*r)->clear_locations();
  }
}

void RateVectorSet::check_duplicate_locations() {
  std::set<bpos> locs = {};
  std::set<bpos> duplicates = {};
  std::cout << "Checking for duplicate locations." << std::endl;
  for(auto r = c.begin(); r != c.end(); ++r) {
    std::cout << (*r)->locations.size() << std::endl;
    for(auto it = (*r)->locations.begin(); it != (*r)->locations.end(); ++it) {
      if(locs.find(*it) == locs.end()) {
	locs.insert(*it);
      } else {
	duplicates.insert(*it);
      }
    }
  }
  std::set<BranchSegment*> b = {};
  std::set<int> p = {};
  std::cout << "Number of non duplicates: " << locs.size() << std::endl;
  for(auto it = locs.begin(); it != locs.end(); ++it) {
    b.insert(it->branch);
    p.insert(it->pos);
  }
  std::cout << "b size: " << b.size() << std::endl;
  std::cout << "p size: " << p.size() << std::endl;
  std::cout << "Number of duplicates: " << duplicates.size() << std::endl;
}

void RateVectorSet::print() {
  for(std::vector<RateVector*>::iterator it = c.begin(); it != c.end(); ++it) {
    (*it)->print();
  }
}

void RateVectorSet::saveToFile(int gen, double l) {
  static int i = -1;
  ++i;
  for(auto it = c.begin(); it != c.end(); ++it) {
    out_file << i << "," << gen << "," << l << "," << (*it)->name << "," << (*it)->state;
    for(auto jt = (*it)->rates.begin(); jt != (*it)->rates.end(); ++jt) {
      out_file << "," << (*jt)->getValue();
    }
    out_file << std::endl;
  }
}

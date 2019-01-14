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
  logLikelihoods = std::vector<double>(20, 0);

  for(int i = 0; i < rates.size(); i++) {
    valueID_to_state[rates[i]->get_ID()] = i;
  }
}

void RateVector::add_location(int pos, BranchSegment* bs) {
  bpos loc = {bs, pos};
  locations.insert(loc);
}

void RateVector::remove_location(int pos, BranchSegment* bs) {
  bpos loc = {bs, pos};
  locations.erase(loc);
}

void RateVector::clear_locations() {
  locations = {};
}

std::unordered_set<bpos> RateVector::get_locations() {
  return(locations);
}

void RateVector::update() {
  /*
   * Updates the substitution counts associtated with each RateParameter(AbstractValue) and the
   * LogLikelihood scores.
   */
  update_counts();
  update_logLikelihoods();
}

void RateVector::update_logLikelihoods() {
  for(int i = 0; i < rates.size(); i++) {
    logLikelihoods[i] = counts[i] * log(rates[i]->getValue());
  }
}

void RateVector::update_single_logLikelihood(int i) {
  int k = valueID_to_state[i];
  logLikelihoods[k] = counts[k] * log(rates[k]->getValue());
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
    logL += logLikelihoods[i];
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

void RateVectorSet::add(RateVector* rv) {
  for(int i = 0; i < rv->rates.size(); i++) {
    // Set up the AbstractValues themselves.
    rv->rates[i]->add_host_vector(rv);
  }
  // Add RateVector to collection in RateVectorSet. 
  c.push_back(rv);
}

void RateVectorSet::get_counts() {
  // Calls get_counts in each rate vector.
  for(auto r = c.begin(); r != c.end(); ++r) {
    (*r)->update();
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

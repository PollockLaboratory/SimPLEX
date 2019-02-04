#include <iostream>

#include "Environment.h"
#include "IO.h"
#include "SubstitutionModel.h"

extern Environment env;
extern IO::Files files;

SubstitutionModel::SubstitutionModel() {
	substitution_model_out = 0;
}

RateVector* SubstitutionModel::selectRateVector(int state) {
	/*
	 * This is a simple function right now but it will become hugely complex.
	 * Given infomation about a BranchSegment and state of interest will return the corresponding rate vector.
	 */
	return(rateVectors[state]);
}

bool SubstitutionModel::SampleParameters() {
	/*
	 * Samples a single parameter within the parameter set();
	 */
	bool sampleType = components.sample();
	return(sampleType);
}

void SubstitutionModel::accept() {
	/*
	 * Accepts the newly sampled parameter set.
	 */
	components.accept();
}

void SubstitutionModel::reject() {
	/*
	 * Rejects the newly sampled parameter set, and undoes the changes from the previous sampling.
	 */
	components.reject();
}

void SubstitutionModel::printParameters() {
	components.print();
}

int SubstitutionModel::getNumberOfParameters() {
	/*
	 * Finds the number of sampleable parameters aka the length of the size of the parameter set.
	 */
	return(components.size());
}

void SubstitutionModel::get_counts() {
  rateVectors.get_counts();  
}

std::vector<RateVector*> SubstitutionModel::get_RateVectors() {
  return(rateVectors.col);
}

void SubstitutionModel::get_current_parameters(std::list<std::pair<RateVector*, int>>& vector_changes) {
  std::list<AbstractComponent*> cur_params = components.get_current_parameters();
  AbstractValue* v;
  for(auto it = cur_params.begin(); it != cur_params.end(); ++it) {
    v = dynamic_cast<AbstractValue*>(*it);
    if(v != NULL) {
      vector_changes.splice(vector_changes.end(), v->get_host_vectors());
    }
  }
}

void SubstitutionModel::saveToFile(int gen, double l) {
  components.saveToFile(gen, l);
  rateVectors.saveToFile(gen, l);
}

void SubstitutionModel::Terminate() {
  delete substitution_model_out;
}

void SubstitutionModel::add_rate_vector(RateVector* v) {
  components.add_rate_vector(v);
  rateVectors.add(v);
}

void SubstitutionModel::finalize() {
  components.Initialize();
  rateVectors.Initialize();
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
  cq.pop(); // Remove first element of queue.
  if(cq.empty()) {
    return(true);
  } else {
    // Add dependancies for new Abstract component at head of the queue.
    for(auto it = sub_model.components.value_to_dependents[cq.front()].begin(); it != sub_model.components.value_to_dependents[cq.front()].end(); ++it) {
      cq.push(*it);
    }
    location = cq.front()->host_vectors.begin();
    location_iter_end = cq.front()->host_vectors.end();
  }
}

SubstitutionModel::iterator::iterator(SubstitutionModel& s, bool e) : sub_model(s), endQ(e) {
  cq = {};
  cq.push(sub_model.components.get_current_parameter());

  // Add dependancies to queue.
  for(auto it = sub_model.components.value_to_dependents[cq.front()].begin(); it != sub_model.components.value_to_dependents[cq.front()].end(); ++it) {
    cq.push(*it);
  }

  location = cq.front()->host_vectors.begin();
  location_iter_end = cq.front()->host_vectors.end();
 
  while(location == location_iter_end) {
    endQ = step_to_next_component();
    if(endQ) {
      return;
    }
  }
}

SubstitutionModel::iterator& SubstitutionModel::iterator::operator++() {
  step_to_next_location();
  return(*this);
}

SubstitutionModel::iterator& SubstitutionModel::iterator::operator++(int i) {
  return(++(*this));
}

const std::pair<RateVector*, int>& SubstitutionModel::iterator::operator*() const {
  return(*location);
}

bool SubstitutionModel::iterator::at_end() const {
  return(endQ);
}

SubstitutionModel::iterator SubstitutionModel::changed_vectors_begin() {
  return(SubstitutionModel::iterator(*this, false));
}

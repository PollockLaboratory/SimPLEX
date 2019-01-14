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
	bool sampleType = parameters.sample();
	return(sampleType);
}

void SubstitutionModel::accept() {
	/*
	 * Accepts the newly sampled parameter set.
	 */
	parameters.accept();
}

void SubstitutionModel::reject() {
	/*
	 * Rejects the newly sampled parameter set, and undoes the changes from the previous sampling.
	 */
	parameters.reject();
}

void SubstitutionModel::printParameters() {
	parameters.print();
}

int SubstitutionModel::getNumberOfParameters() {
	/*
	 * Finds the number of sampleable parameters aka the length of the size of the parameter set.
	 */
	return(parameters.size());
}

void SubstitutionModel::get_counts() {
  rateVectors.get_counts();  
}


double SubstitutionModel::get_substitution_logLikelihood() {
  double logL = 0.0;
  for(auto it = rateVectors.c.begin(); it != rateVectors.c.end(); ++it) {
    //(*it)->update_logLikelihoods();
    logL += (*it)->get_logLikelihood();
  }
  return(logL);
}

void SubstitutionModel::clear_locations() {
  rateVectors.clear_locations();
}

void SubstitutionModel::check_duplicate_locations() {
  rateVectors.check_duplicate_locations();
}

std::list<AbstractValue*> SubstitutionModel::get_current_parameters() {
	return(parameters.get_current_parameters());
}

void SubstitutionModel::saveToFile(int gen, double l) {
  parameters.saveToFile(gen, l);
  rateVectors.saveToFile(gen, l);
}

void SubstitutionModel::Terminate() {
  delete substitution_model_out;
}

void SubstitutionModel::add_rate_vector(RateVector* v) {
  parameters.add_rate_vector(v);
  rateVectors.add(v);
}

void SubstitutionModel::finalize() {
  parameters.Initialize();
  rateVectors.Initialize();
}

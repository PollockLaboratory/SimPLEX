#include <iostream>
#include <cmath>
#include <sys/times.h>
#include <chrono>

#include "Model.h"
#include "Trees/TreeParser.h"
#include "SubstitutionModels/SubstitutionModel.h"

#include "Environment.h"
#include "IO.h"

extern Environment env;
extern IO::Files files;

/*  The tree needs the names_to_sequence map. The substitution model might need the empirical frequencies. */

// using namespace std;

Model::Model() {
  /*
   * The default contructor.
   */
  tree = NULL;
  ready = true;
}

void Model::Initialize(IO::RawTreeNode* &raw_tree, SequenceAlignment* &MSA, SubstitutionModel* &sm) {
  /*
   * Initialize the model class.
   * There are two main components within the model class:
   * - the tree class - contains the tree topology as well as the sequences.
   * - the substitution model class - which contains all the rate matrices.
   */
  substitution_model = sm;
  num_parameters = substitution_model->getNumberOfParameters();

  tree = new Tree();
  tree->Initialize(raw_tree, MSA, substitution_model);

  // Sort out counts.
  counts = SubstitutionCounts(substitution_model->get_RateVectors(), tree->get_branch_lengths());
  tree->update_counts(counts);
  counts.print();
}

// Sampling
bool Model::SampleTree() {
  if(ready) {
    bool t = tree->sample();
    counts = SubstitutionCounts(substitution_model->get_RateVectors(), tree->get_branch_lengths()); // Empty list of counts
    tree->update_counts(counts);

    return(t);
  } else {
    std::cout << "Error: Attempt to sample tree before accepting previous changes." << std::endl;
    exit(EXIT_FAILURE);
  }
}

bool Model::SampleSubstitutionModel() {
  if(ready) {
    return(substitution_model->SampleParameters());
    ready = false;
  } else {
    std::cout << "Error: Attempt to sample parameter before accepting previous changes." << std::endl;
    exit(EXIT_FAILURE);
  }
}

void Model::accept() {
  substitution_model->accept();
  ready = true;
}

void Model::reject() {
  substitution_model->reject();
  logL -= delta_logL;
  ready = true;
}

// Likelihood Calculations.
double Model::CalculateLikelihood() {
  /*
   * Calculates the likelihood of the current tree and substitution model.
   */

  logL_waiting = 0.0;
  logL_subs = 0.0;
  logL = 0.0;
  float t;
  int num0subs;
  int num1subs;

  double u = substitution_model->get_u();

  for(auto it = counts.subs_by_branch.begin(); it != counts.subs_by_branch.end(); ++it) {
    t = it->first;
    num0subs = it->second.num0subs;
    num1subs = it->second.num1subs;
    logL_waiting += num0subs * log(1/(1 + u*t)) + num1subs * log(t/(1 + u*t));
  }

  for(auto it = counts.subs_by_rateVector.begin(); it != counts.subs_by_rateVector.end(); ++it) {
    RateVector* rv = it->first;
    std::vector<int> C_xy = it->second;
    for(int i = 0; i < env.num_states; i++) {
      logL_subs += C_xy[i] * log((*rv)[i]);
    }
  }

  logL = logL_waiting + logL_subs;

  return(logL);
}

double Model::updateLikelihood(){
  /*
   * Rather recalculating the full likelihood, will modify logLikelihood only by what has changed. 
   */
  delta_logL = 0.0;
 
  RateVector* rv;
  int C_xy;

  for(auto it = substitution_model->changed_vectors_begin(); it.at_end() == false; ++it) {
    rv = (*it).rv;
    C_xy = counts.subs_by_rateVector[rv][(*it).pos];
    delta_logL += C_xy * log(rv->get_rate_ratio((*it).pos)); // Should be ratio;
  }

  logL += delta_logL;
  return(logL);
}

// Printing/Recording
void Model::RecordState(int gen, double l) {
	/*
	 * Records the state of both the tree and the substitution model.
	 */
	tree->record_state(gen, l);
	substitution_model->saveToFile(gen, l);
}

void Model::print() {

}

void Model::printParameters() {
  substitution_model->printParameters();
  //tree->printCounts();
}

// Tidying up
void Model::Terminate() {
	/*
	 * Just terminates substitution_model.
	 */
  // Save the tree data.
  tree->record_tree();
  substitution_model->Terminate();
}


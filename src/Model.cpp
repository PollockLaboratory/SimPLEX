#include <iostream>
#include <cmath>
#include <sys/times.h>
#include <chrono>

#include "Model.h"
#include "Trees/TreeParser.h"
#include "SubstitutionModels/SubstitutionModelTypes.h"
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
  SubstitutionModel* substitution_model = NULL;
  ready = true;
}

Model::~Model() {
  /*
   * The destructor function.
   * Note: we should not be destructing manually.
   */
  delete tree;
  delete substitution_model;
}

SubstitutionModel* Model::InitializeSubstitutionModel(int num_sites, vector<string> states) {
  /*
   * This is in the substituionModelType.h, can we make this more explicit?
   */
  std::cout << "Initializing Substitution Model: ";
  SubstitutionModel* substitution_model = GetSubstitutionModel(); // In SubstituionModelTypes.h
  substitution_model->Initialize(num_sites, states);
  std::cout << std::endl;
  return(substitution_model);
}

void Model::Initialize(IO::RawTreeNode* &raw_tree, SequenceAlignment* &MSA) {
  /*
   * Initialize the model class.
   * There are two main components within the model class:
   * - the tree class - contains the tree topology as well as the sequences.
   * - the substitution model class - which contains all the rate matrices.
   */
  u = env.u;
  int num_sites = (MSA->taxa_names_to_sequences).begin()->second.size();
  substitution_model = InitializeSubstitutionModel(num_sites, MSA->states);
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

  for(auto it = counts.subs_by_branch.begin(); it != counts.subs_by_branch.end(); ++it) {
    t = it->first;
    num0subs = it->second.first;
    num1subs = it->second.second;
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

  // double l = tree->calculate_likelihood();
  return(logL);
}

double Model::updateLikelihood(){
  /*
   * Rather recalculating the full likelihood, will modify logLikelihood only by what has changed. 
   */
  delta_logL = 0.0;

  int time_taken;
  std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
  
  std::list<std::pair<RateVector*, int>> vector_changes = substitution_model->get_current_parameters();

  time_taken = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - start).count();
  std::cout << "Time: " << time_taken << std::endl;

  for(auto it = vector_changes.begin(); it != vector_changes.end(); ++it) {
    RateVector* rv = it->first;
    std::vector<int> C_xy = counts.subs_by_rateVector[rv];
    delta_logL += C_xy[it->second] * log(rv->get_rate_ratio(it->second)); // Should be ratio;
  }

  logL += delta_logL;
  // double l = tree->update_likelihood();
  // std::cout << "L1: " << l << " L2: " << logL << std::endl;
  return(logL);
}

// double Model::PartialCalculateLikelihood(const double lnL) {
// return(tree->partial_calculate_likelihood());
//}

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


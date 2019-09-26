#include <iostream>
#include <cmath>
#include <sys/times.h>
#include <chrono>

#include "Model.h"
#include "ModelParts/SubstitutionModels/SubstitutionModel.h"

#include "Environment.h"
#include "IO/Files.h"
#include "IO/TreeParser.h"

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

void Model::Initialize(IO::RawTreeNode* &raw_tree, IO::RawMSA* &raw_msa, IO::raw_substitution_model* &raw_sm) {

  UniformizationConstant* u = new UniformizationConstant();
  components.add_parameter(u);
  
  // Substitution model.
  substitution_model = new SubstitutionModel(u);
  substitution_model->from_raw_model(raw_sm);

  std::list<AbstractComponent*> sm_parameters = substitution_model->get_all_parameters();
  for(auto it = sm_parameters.begin(); it != sm_parameters.end(); ++it) {
    components.add_parameter(*it);
  }

  components.refresh_dependencies();
  u->set_initial();
  components.Initialize();

  // Sequences.
  const States* states = substitution_model->get_states();

  SequenceAlignment* MSA = new SequenceAlignment(states);
  MSA->Initialize(raw_msa);

  // Tree.
  tree = new Tree();
  tree->Initialize(raw_tree, MSA, substitution_model);

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
    std::cerr << "Error: Attempt to sample tree before accepting previous changes." << std::endl;
    exit(EXIT_FAILURE);
  }
}

bool Model::SampleSubstitutionModel() {
  if(ready) {
    return(components.sample());
    ready = false;
  } else {
    std::cerr << "Error: Attempt to sample parameter before accepting previous changes." << std::endl;
    exit(EXIT_FAILURE);
  }
}

void Model::accept() {
  components.accept();
  ready = true;
}

void Model::reject() {
  components.reject();
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

  for(auto it = substitution_model->modified_begin(components.get_current_parameter()); it.at_end() == false; ++it) {
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
	components.saveToFile(gen, l);
	substitution_model->saveToFile(gen, l);
}

void Model::print() {
}

void Model::printParameters() {
  components.print();
  //substitution_model->printParameters();
  //tree->printCounts();
}

// Tidying up
void Model::Terminate() {
  /*
   * Just terminates substitution_model.
   */
  // Save the tree data.
  tree->record_tree();
}


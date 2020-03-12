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

Model::Model() {
  /*
   * The default contructor.
   */
  ready = true;
}

Valuable* Model::create_uniformization_constant() {
  Valuable* val =  nullptr;
  if(env.get<bool>("UNIFORMIZATION.dynamic")) {
    UniformizationConstant* u = new UniformizationConstant(env.get<double>("UNIFORMIZATION.initial_value"));
    //u->set_initial();
    components.add_parameter(u, env.get<int>("UNIFORMIZATION.sample_frequency"));
    val = u; 
  } else {
    FixedFloat* f = new FixedFloat("U", env.get<double>("UNIFORMIZATION.initial_value"));
    components.add_parameter(f);
    val = f;
  }
  return(val);
}

void Model::Initialize(IO::RawTreeNode* &raw_tree, IO::RawMSA* &raw_msa, IO::raw_substitution_model* &raw_sm) {

  std::cout << "Initializing model:" << std::endl;

  Valuable* u = create_uniformization_constant();

  // Substitution model.
  std::cout << "\tConstructing Substitution model." << std::endl;
  substitution_model = new SubstitutionModel(u);
  substitution_model->from_raw_model(raw_sm);

  // Add all substitution parameters to component set.
  std::list<AbstractComponent*> sm_parameters = substitution_model->get_all_parameters();
  for(auto it = sm_parameters.begin(); it != sm_parameters.end(); ++it) {
    components.add_parameter(*it);
  }

  // Tree.
  std::cout << "\tConstructing tree." << std::endl;
  AncestralStatesParameter* tp = new AncestralStatesParameter();
  tp->Initialize(raw_tree, raw_msa, substitution_model);

  std::cout << "\tAdding Uniformization constant." << std::endl;
  UniformizationConstant* u_param = dynamic_cast<UniformizationConstant*>(u);
  if(u_param and env.get<bool>("UNIFORMIZATION.refresh_tree_on_update")) {
      tp->add_dependancy(u_param);
  }

  components.add_parameter(tp, env.get<int>("MCMC.tree_sample_frequency"));

  
  std::cout << "\tPreparing substitution counts." << std::endl;
  CountsParameter* cp = new CountsParameter(&counts, tp);
  components.add_parameter(cp);

  components.Initialize();

  std::cout << "\tSetting initial parameter states." << std::endl;
  // Set initial tree state.
  tp->sample();
  components.reset_dependencies();

  counts.print();
  std::cout << "Model succesfully constructed." << std::endl;
}

// Sampling
sample_status Model::sample() {
  if(ready) {
    //std::cout << std::endl << "New sample." << std::endl;
    return(components.sample());
    ready = false;
  } else {
    std::cerr << "Error: Attempt to sample next parameter before accepting previous changes." << std::endl;
    exit(EXIT_FAILURE);
  }
}

void Model::accept() {
  components.accept();
  ready = true;
}

void Model::reject() {
  components.reject();
  //logL = previous_logL;
  ready = true;
}

// Likelihood Calculations.
double Model::CalculateLikelihood() {
  /*
   * Calculates the likelihood of the current tree and substitution model.
   */
  
  double logL_waiting = 0.0;
  double logL_subs = 0.0;

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
  
  return(logL_waiting + logL_subs);
}

double Model::CalculateChangeInLikelihood(){
  /*
   * Rather recalculating the full likelihood, will modify logLikelihood only by what has changed. 
   */

  double delta_logL = 0.0;
 
  RateVector* rv;
  int C_xy;

  //std::cout << std::endl;
  for(auto it = substitution_model->modified_begin(components.get_current_parameter());
      it.at_end() == false; ++it) { 
    rv = (*it).rv;
    C_xy = counts.subs_by_rateVector[rv][(*it).pos];
    delta_logL += C_xy * log(rv->get_rate_ratio((*it).pos));
    //std::cout << "[" << rv->get_rate_ratio((*it).pos) << " " << C_xy * log(rv->get_rate_ratio((*it).pos)) << "] ";
    //std::cout << "RV: " << rv->get_name() << " count: " << C_xy << " ratio: " << rv->get_rate_ratio((*it).pos) << " log: " << log(rv->get_rate_ratio((*it).pos)) << " " << delta_logL << std::endl;
  }

  //std::cout << "Delta: " << delta_logL << std::endl;
  return(delta_logL);
}

// Printing/Recording
void Model::RecordState(int gen, double l) {
	/*
	 * Records the state of both the tree and the substitution model.
	 */
	components.saveToFile(gen, l);
	substitution_model->saveToFile(gen, l);
}

void Model::print() {
}

void Model::printParameters() {
  components.print();
  //substitution_model->printParameters();
}

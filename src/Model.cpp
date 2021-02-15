#include <iostream>
#include <map>

#include "Model.h"

#include "Environment.h"
#include "IO/Files.h"

#include "ModelParts/SubstitutionModels/Parameters.h"

extern Environment env;
extern IO::Files files;

/*  The tree needs the names_to_sequence map. The substitution model might need the empirical frequencies. */

Model::Model() {
  /*
   * The default contructor.
   */
  ready = true;
  components.set_counts(&counts);
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
    f->hide();
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

  // Primary States.
  const States* states = substitution_model->get_states("primary");
  std::cout << "Primary states: ";
  print_States(*states);

  SequenceAlignment* MSA = new SequenceAlignment("primary",
						 env.get<std::string>("OUTPUT.sequences_out_file"),
						 env.get<std::string>("OUTPUT.substitutions_out_file"),
						 states);
  MSA->Initialize(raw_msa);

  std::map<std::string, SequenceAlignment*> hidden_states_MSAs = {};
  // Secondary/Hidden states.
  std::map<std::string, std::list<std::string>> all_states = raw_sm->get_all_states();
  for(auto it = all_states.begin(); it != all_states.end(); ++it) {
    if(it->first != "primary") {
      const States *secondary_states = substitution_model->get_states(it->first);
      std::cout << "Secondary states: ";
      print_States(*secondary_states);
      SequenceAlignment* secondary_MSA = new SequenceAlignment(it->first,
							       raw_sm->states_seqs_output_files[it->first],
							       raw_sm->states_subs_output_files[it->first],
							       secondary_states);
      secondary_MSA->Initialize(raw_sm->get_hidden_states_data(it->first));

      hidden_states_MSAs[it->first] = secondary_MSA;

      if(MSA->match_structure(secondary_MSA) == false) {
	std::cerr << "Error: sequence alignment for domain \"" << it->first << "\" does not have a matching structure to the primary alignment." << std::endl;
	exit(EXIT_FAILURE);
      }
    }
  }

  // Set environment variables.
  env.num_states = states->n;
  env.state_to_integer = states->state_to_int;
  env.integer_to_state = states->int_to_state;

  // Tree.
  std::cout << "\tConstructing tree." << std::endl;

  tree = new Tree();
  tree->Initialize(raw_tree);
  tree->connect_substitution_model(substitution_model);

  tree->MSA = MSA;
  tree->configureBranches(tree->root, MSA->n_columns, all_states);

  // Configuring sequences.
  MSA->syncWithTree(tree);

  RateVectorAssignmentParameter* rvap = new RateVectorAssignmentParameter(tree);

  unsigned int ctr = 0;
  for(auto it = hidden_states_MSAs.begin(); it != hidden_states_MSAs.end(); ++it) {
    it->second->syncHiddenWithTree(it->first, ctr, tree);
    ctr++;
    SequenceAlignmentParameter* hidden_parameter = new SequenceAlignmentParameter(it->second);
    components.add_parameter(hidden_parameter, env.get<int>("MCMC.hidden_sample_frequency"));
    rvap->add_dependancy(hidden_parameter);
    msa_parameters.push_back(hidden_parameter);
  }

  SequenceAlignmentParameter* sp = new SequenceAlignmentParameter(MSA);
  msa_parameters.push_back(sp);

  components.add_parameter(sp, env.get<int>("MCMC.tree_sample_frequency"));

  // This should be closer to RateVector.Initialize() call.
  substitution_model->organizeRateVectors(MSA->numCols());

  rvap->add_dependancy(sp);

  components.add_parameter(rvap);

  std::cout << "\tAdding Uniformization constant." << std::endl;
  UniformizationConstant* u_param = dynamic_cast<UniformizationConstant*>(u);
  if(u_param and env.get<bool>("UNIFORMIZATION.refresh_tree_on_update")) {
      sp->add_dependancy(u_param);
  }
 
  std::cout << "\tPreparing substitution counts." << std::endl;
  cp = new CountsParameter(&counts, tree, all_states);
  cp->add_dependancy(rvap);
  components.add_parameter(cp);

  components.Initialize();

  std::cout << "\tSetting initial parameter states." << std::endl;
  // Set initial tree state.
  sp->sample();
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

  //std::cout << "Counts: ";
  for(auto it = counts.subs_by_rateVector.begin(); it != counts.subs_by_rateVector.end(); ++it) {
    RateVector* rv = it->first;
    //std::cout << rv->get_name() << " [ ";
    std::vector<int> C_xy = it->second;
    for(int i = 0; i < rv->rates.size(); i++) {
      // std::cout << C_xy[i] << "|" << log((*rv)[i]) << " ";
      logL_subs += C_xy[i] * log((*rv)[i]);
    }
    //std::cout << "]" << std::endl;
  }

  //std::cout << logL_waiting << " " << logL_subs << std::endl;
 
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
	components.save_to_file(gen, l);
	substitution_model->saveToFile(gen, l, counts.subs_by_rateVector);
	for(auto it = msa_parameters.begin(); it != msa_parameters.end(); ++it) {
	  (*it)->save_to_file(gen, l);
	}
	// Substitutions should be a call to the MSA.
	//tree->record_substitutions(gen, l);
	cp->save_to_file(gen, l); 
}

void Model::print() {
}

void Model::printParameters() {
  components.print();
  //substitution_model->printParameters();
}

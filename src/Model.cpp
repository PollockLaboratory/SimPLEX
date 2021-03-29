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

void Model::Initialize(IO::RawTreeNode* &raw_tree, IO::raw_substitution_model* &raw_sm) {

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

  std::list<std::string> tip_sequences = getRawTreeNodeTipNames(raw_tree);

  // All states.
  unsigned int n_cols = 0;
  std::map<std::string, SequenceAlignment*> all_MSAs = {};
  std::map<std::string, std::list<std::string>> all_states = raw_sm->get_all_states();
  for(auto it = all_states.begin(); it != all_states.end(); ++it) {
      const States *secondary_states = substitution_model->get_states(it->first);
      std::cout << "Reading new states domain: ";
      print_States(*secondary_states);
      SequenceAlignment* secondary_MSA = new SequenceAlignment(it->first,
							       raw_sm->states_seqs_output_files[it->first],
							       raw_sm->states_subs_output_files[it->first],
							       secondary_states);
      secondary_MSA->Initialize(raw_sm->get_state_data(it->first));

      // Validate MSAs.
      // Check node names are consistant between tree and sequence alignment.
      // TODO - Check pattern of gaps is consistant.
      if(secondary_MSA->validate(tip_sequences, all_MSAs) != true) {
	std::cerr << "Error: sequence alignment for domain \"" << it->first << "\" cannot be validated." << std::endl;
	exit(EXIT_FAILURE);
      }

      all_MSAs[it->first] = secondary_MSA;

      n_cols = secondary_MSA->numCols();
  }

  // COnfiguring the Tree.
  std::cout << "\tConstructing tree." << std::endl;

  tree = new Tree();
  tree->Initialize(raw_tree);
  tree->connect_substitution_model(substitution_model);
  tree->configureBranches(tree->root, n_cols, all_states);

  // Configuring sequences.
  RateVectorAssignmentParameter* rvap = new RateVectorAssignmentParameter(tree);

  UniformizationConstant* u_param = dynamic_cast<UniformizationConstant*>(u);
  unsigned int ctr = 0;
  for(auto it = all_MSAs.begin(); it != all_MSAs.end(); ++it) {
    it->second->syncWithTree(it->first, ctr, tree);
    ctr++;
    SequenceAlignmentParameter* msa_parameter = new SequenceAlignmentParameter(it->second);
    if(u_param and env.get<bool>("UNIFORMIZATION.refresh_tree_on_update")) {
      msa_parameter->add_dependancy(u_param);
    }
    components.add_state_parameter(msa_parameter, env.get<int>("MCMC.alignment_sample_frequency"));
    rvap->add_dependancy(msa_parameter);
    msa_parameters.push_back(msa_parameter);
  }

  // This should be closer to RateVector.Initialize() call.
  substitution_model->organizeRateVectors();

  components.add_parameter(rvap);

  std::cout << "\tPreparing substitution counts." << std::endl;
  cp = new CountsParameter(&counts, tree, all_states);
  cp->add_dependancy(rvap);
  components.add_parameter(cp);

  components.Initialize();

  std::cout << "\tSetting initial parameter states." << std::endl;
  // Set parameter states.
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
    for(unsigned int i = 0; i < rv->rates.size(); i++) {
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

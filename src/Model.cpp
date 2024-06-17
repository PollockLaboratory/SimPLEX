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
  /*
   * This create the uniformization constant as a model parameter.
   * This is not strictly necessary as it is a static value, however this was more complicated
   * when there was a possibility that the uniformization constant could change.
   * This would be a cool feature, but it appears more trouble than its worth at this point.
   */
  FixedFloat* f = new FixedFloat("U", env.get<double>("UNIFORMIZATION.initial_value"));
  f->hide();
  components.add_parameter(f);
  return(f);
}

SequenceAlignment* initialize_dynamic_alignment(std::string name, const IO::RawMSA* raw_msa, const States* states, std::string seqs_out_file, std::string subs_out_file) {
  SequenceAlignment* MSA = new SequenceAlignment(name, seqs_out_file, subs_out_file, states); 
  MSA->initialize_dynamic(*raw_msa);

  return MSA;
}

SequenceAlignment* initialize_site_static_alignment(std::string name, const IO::RawMSA* raw_msa, const States* states, std::string seqs_out_file, std::string subs_out_file) {
  SequenceAlignment* MSA = new SequenceAlignment(name, seqs_out_file, subs_out_file, states); 
  MSA->initialize_site_static(*raw_msa);

  return MSA;
}

std::optional<unsigned int> find_column_count(const std::map<std::string, SequenceAlignment*>& all_MSAs) {
  std::optional<unsigned int> n_cols;
  for(auto it = all_MSAs.begin(); it != all_MSAs.end(); ++it) {
    if (n_cols) {
      // If the column number is not consistant return null option.
      if (n_cols.value() != it->second->n_cols()) return std::nullopt;
    } else {
      n_cols = std::optional<unsigned int>(it->second->n_cols());
    }
  }
  return n_cols;
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

  // STATE DOMAINS / DATA 
  std::map<std::string, SequenceAlignment*> all_MSAs = {};
  for (const auto& [state_domain, states] : substitution_model->get_all_states()) {
    std::cout << "\t\tReading new states domain \'" << state_domain << "\': "; print_States(states);

    IO::StateData state_data = raw_sm->get_state_data(state_domain);
    switch(state_data.tag) {
    case IO::StateData::Tag::DYNAMIC: {
      const IO::RawMSA* raw_msa = state_data.data.dynamic;
      SequenceAlignment* MSA = initialize_dynamic_alignment(state_domain, raw_msa, &states, raw_sm->states_seqs_output_files[state_domain], raw_sm->states_subs_output_files[state_domain]);
 
      // Validate MSAs.
      // Check node names are consistant between tree and sequence alignment.
      // TODO - Check pattern of gaps is consistant.
      if(MSA->validate(tip_sequences, all_MSAs) != true) {
        std::cerr << "Error: sequence alignment for domain \"" << state_domain << "\" cannot be validated." << std::endl;
        exit(EXIT_FAILURE);
      }

      all_MSAs[state_domain] = MSA;
      break;
    }
    case IO::StateData::Tag::SITE_STATIC: {
      const IO::RawMSA* raw_msa = state_data.data.site_static;
      SequenceAlignment* MSA = initialize_site_static_alignment(state_domain, raw_msa, &states, raw_sm->states_seqs_output_files[state_domain], raw_sm->states_subs_output_files[state_domain]);

      // VALIDATE
      all_MSAs[state_domain] = MSA;

      break;
    }
    default: {
      std::cout << "Error: state data type not recognised." << std::endl;
      exit(EXIT_FAILURE);
      break;
    }
    }
  }

  std::optional<unsigned int> n_col_opt = find_column_count(all_MSAs);
  if (not n_col_opt.has_value()) {
    std::cout << "Error: MSA column counts are not equal across states." << std::endl;
    exit(EXIT_FAILURE);
  }

  // Configuring the Tree.
  std::cout << "\tConstructing tree." << std::endl;

  std::list<std::string> state_domain_names = {};
  for (const auto& [state_domain, states] : substitution_model->get_all_states()) state_domain_names.push_back(state_domain);
 
  tree = new Tree();
  tree->Initialize(raw_tree);
  tree->connect_substitution_model(substitution_model);
  tree->configure_branches(tree->root, n_col_opt.value(), state_domain_names);

  // Configuring sequences.
  unsigned int n_sample = env.get<unsigned int>("MCMC.position_sample_count");
  RateVectorAssignmentParameter* rvap = new RateVectorAssignmentParameter(tree);

  for(const auto& [state_domain, msa] : all_MSAs) {
    switch (msa->tag) {
    case SequenceAlignment::Tag::DYNAMIC: {
      //std::cout << state_domain << "->DYNAMIC" << std::endl;
      msa->syncWithTree(state_domain, tree);
      SequenceAlignmentParameter* msa_parameter = new SequenceAlignmentParameter(msa, n_sample);
      msa_parameters.push_back(msa_parameter);

      components.add_state_parameter(msa_parameter, env.get<int>("MCMC.alignment_sample_frequency"));
      rvap->add_dependancy(msa_parameter);
      break;
    }
    case SequenceAlignment::Tag::SITE_STATIC: {
      //std::cout << state_domain << "->SITE_STATIC" << std::endl;
      msa->syncWithTree(state_domain, tree);
      this->substitution_model->mark_static_state(state_domain);

      break;
    }
    default: {
      std::cout << "Error: sequence alignment type not implemented." << std::endl;
      exit(EXIT_FAILURE);
      break;
    }
    }
  }

  // This should be closer to RateVector.Initialize() call.
  substitution_model->organizeRateVectors();

  components.add_parameter(rvap);

  std::cout << "\tPreparing substitution counts." << std::endl;

  std::map<std::string, std::list<std::string>> all_state_domains = raw_sm->get_all_states();
  cp = new CountsParameter(&counts, tree, all_state_domains);
  cp->add_dependancy(rvap);
  components.add_parameter(cp);

  components.Initialize();

  std::cout << "\tSetting initial parameter states." << std::endl;

  // Set parameter states.
  components.reset_dependencies();

  std::cout << "Model succesfully constructed." << std::endl << std::endl;

  std::cout << "Initial substitution counts:" << std::endl;
  counts.print();
}

// Sampling
sample_status Model::sample() {
  if(ready) {
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

  int num0subs;
  int num1subs;

  double u = substitution_model->get_u();

  for(const auto& [branch_length, branch_counts] : this->counts.subs_by_branch) {
    num0subs = branch_counts.num0subs;
    num1subs = branch_counts.num1subs;

    logL_waiting += num0subs * log(1/(1 + u*branch_length)) + num1subs * log(branch_length/(1 + u*branch_length));
  }

  for(auto it = counts.subs_by_rateVector.begin(); it != counts.subs_by_rateVector.end(); ++it) {
    RateVector* rv = it->first;
    std::vector<double> C_xy = it->second;
    for(unsigned int i = 0; i < rv->rates.size(); i++) {
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

  for(auto it = substitution_model->modified_begin(components.get_current_parameter());
      it.at_end() == false; ++it) { 
    rv = (*it).rv;
    C_xy = counts.subs_by_rateVector[rv][(*it).pos];
    delta_logL += C_xy * log(rv->get_rate_ratio((*it).pos));
  }

  return(delta_logL);
}

// Printing/Recording
void Model::RecordState(uint128_t gen, double l) {
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

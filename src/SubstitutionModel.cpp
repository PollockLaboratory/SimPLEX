#include "SubstitutionModel.h"

#include <iostream>

#include "Options.h"

extern Options options;

SubstitutionModel::SubstitutionModel() {
	id = 0;
	is_constant = false;
	substitution_model_out = 0;
}

SubstitutionModel::~SubstitutionModel() {

}

void SubstitutionModel::Initialize() {
//	std::cout << "Initializing Substitution Model" << std::endl;

	is_constant = options.constant_substitution_models.front();
	// This allows recycling of the last value
	if (options.constant_substitution_models.size() > 1) options.constant_substitution_models.pop();
}

// Is this the best place for these? Perhaps below?
static const int initialize_by_sampling = 0;
static const int initialize_from_file = 1;
static const int initialize_by_asking_for_file = 2;

/**
 * I like that this is a unified interface. Since it is unified, it should
 * take as much information as required by the most complex model I have.
 * That means passing in number_of_sites, residue_to_integer, and
 * integer_to_residue.
 *
 */
void SubstitutionModel::InitializeState() {
	std::cout << "Initializing model " << IdToString();
	int initialization_type =
			options.substitution_models_initialization_type.front();
	// This allows for recycling of the final model initialization type
	if (options.substitution_models_initialization_type.size() > 1)
		options.substitution_models_initialization_type.pop();

	if (initialization_type == initialize_by_sampling) {
		std::cout << " randomly" << std::endl;
		SampleParameters();
	} else if (initialization_type == initialize_from_file) {
		// What if initializing state from file requires more information than
		// simply the file? What if it requires the residues too?
		// Perhaps the file should give the residues too.
		std::string state_in_file =
				options.substitution_models_initialization_files.front();
		options.substitution_models_initialization_files.pop();

		std::cout << " from file " << state_in_file <<  std::endl;
// I want to open the file here and pass in the ifstream. All the substitution
		// models will have to open the file so I should do it here instead.
		InitializeStateFromFile(state_in_file.c_str());

	} else if (initialization_type == initialize_by_asking_for_file) {
		std::cout << std::endl;
		std::cout << "What is the initialization file for substitution model "
				<< IdToString() << "?" << std::endl;

		std::string state_in_file;
		std::cin >> state_in_file;

		InitializeStateFromFile(state_in_file);
	}
}

//These are purely virtual for now.
//void SubstitutionModel::SampleParameters() {}
//void SubstitutionModel::RecordState() {}

void SubstitutionModel::Terminate() {
	delete substitution_model_out;
}

std::string SubstitutionModel::IdToString() {
	std::ostringstream s;
	s << id;
	return s.str();
}

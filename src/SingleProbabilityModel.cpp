#include "SingleProbabilityModel.h"

#include <iostream>
#include <cstdlib>

#include "Options.h"

extern Options options;

extern double Random();

int SingleProbabilityModel::number_of_single_probability_models = 0;

SingleProbabilityModel::SingleProbabilityModel() {
	id = number_of_single_probability_models;
	number_of_single_probability_models++;

	substitution_probability = 0;
}

SingleProbabilityModel::SingleProbabilityModel(
		const SingleProbabilityModel& single_probability_model) :
		SubstitutionModel(single_probability_model) {

	substitution_probability =
			single_probability_model.substitution_probability;
}

SingleProbabilityModel::~SingleProbabilityModel() {
//	std::cout << "Single probability model deconstructor" << std::endl;
//	delete substitution_model_out;
}

SingleProbabilityModel* SingleProbabilityModel::Clone() {
//	std::cout << "cloning single probability model" << std::endl;
	return new SingleProbabilityModel(*this);
}

void SingleProbabilityModel::Initialize(int number_of_sites,
		std::vector<std::string> states) {
//	std::cout << "Initializing Single Probability Model" << std::endl;
	SubstitutionModel::Initialize(); // sets is_constant

	InitializeState(); // Can call InitializeStateFromFile()
	InitializeOutputStream();
	if (is_constant) {
		is_constant = false;
		RecordState();
		is_constant = true;
	}
}

void SingleProbabilityModel::InitializeStateFromFile(
		std::string state_in_file) {
	std::ifstream state_in(state_in_file.c_str());

	std::string header;
	state_in >> header;
	if (header != "Substitution_Probability") {
		std::cerr
				<< "Check header on substitution probability model initialization file for model "
				<< IdToString() << std::endl;
		std::cerr << "The header should be \"Substitution_Probability\""
				<< std::endl;
		std::exit(-1);
	}

	state_in >> substitution_probability;
}

void SingleProbabilityModel::InitializeOutputStream() {
	std::string substitution_model_out_file = "single_probability_model_"
			+ IdToString();
	options.PrependOutputDirectory(substitution_model_out_file);
	substitution_model_out = new std::ofstream(
			substitution_model_out_file.c_str());
	*substitution_model_out << "Substitution_Probability" << std::endl;
}

/**
 * Notice that this is independent of the current substitution_probability.
 * This is the simplest (non-constant) way to sample a probability.
 * This samples from the uniform prior and not from the posterior or
 * marginal.
 *
 */
void SingleProbabilityModel::SampleParameters() {
//	std::cout << "Sampling single probability model parameters" << std::endl;
	if (not is_constant) {
		// This is uniform prior Metropolis-Hastings sampling
		// It also forgets the previous state. Is that preferrable?
		substitution_probability = Random();
	}
}

/**
 * Only record the state if the model is not constant. If it is constant,
 * then the parameters will not change and there is no reason to print it
 * over and over.
 *
 */
void SingleProbabilityModel::RecordState() {
	if (not is_constant)
		*substitution_model_out << substitution_probability << std::endl;
}

double SingleProbabilityModel::SubstitutionProbability(int ancestral_state,
		int descendent_state, int site, double branch_length) {
	double probability = 0;

	if (ancestral_state != descendent_state)
		probability = substitution_probability;
	else
		probability = 1 - substitution_probability;

	return probability;
}

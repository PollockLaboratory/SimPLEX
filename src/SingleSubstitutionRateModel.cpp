#include "SingleSubstitutionRateModel.h"

#include <cstdlib> // For exit()
#include <iostream>
#include <cmath> // For exp()

#include "Options.h"
extern Options options;

extern double Random();

int SingleSubstitutionRateModel::number_of_single_substitution_rate_models = 0;

SingleSubstitutionRateModel::SingleSubstitutionRateModel() {
	id = number_of_single_substitution_rate_models;
	number_of_single_substitution_rate_models++;
	substitution_rate = 0;
}

SingleSubstitutionRateModel::~SingleSubstitutionRateModel() {
//	delete substitution_model_out;
}

SingleSubstitutionRateModel* SingleSubstitutionRateModel::Clone() {
//	std::cout << "cloning single substitution rate model" << std::endl;
	return new SingleSubstitutionRateModel(*this);
}

void SingleSubstitutionRateModel::InitializeOutputStream() {
	std::string substitution_model_out_file = "single_substitution_rate_model_"
			+ IdToString();
	options.PrependOutputDirectory(substitution_model_out_file);
	substitution_model_out = new std::ofstream(
			substitution_model_out_file.c_str());
	*substitution_model_out << "Substitution_Rate" << std::endl;
}

void SingleSubstitutionRateModel::Initialize(int number_of_sites,
		std::vector<std::string> states) {
	std::cout << "Initializing Single Substitution Rate Model" << std::endl;
	SubstitutionModel::Initialize(); // sets is_	constant

	InitializeState(); // Can call SampleParameters or InitializeStateFromFunction

	InitializeOutputStream();

}

void SingleSubstitutionRateModel::InitializeStateFromFile(
		std::string state_in_file) {
	std::ifstream state_in(state_in_file.c_str());

	std::string header;
	state_in >> header;
	if (header != "Substitution_Rate") {
		std::cerr
				<< "Check header on single substitution rate model initialization file for model "
				<< IdToString() << std::endl;
		std::cerr << "The header should be \"Substitution_Rate\"" << std::endl;
		exit(-1);
	}

	state_in >> substitution_rate;
}

/**
 * Notice that this is independent of the current substitution_rate.
 * This is the simplest (non-constant) way to sample a variable.
 * This samples from the uniform prior and not from the posterior or
 * marginal.
 *
 */
static double step_size = 0.1;
void SingleSubstitutionRateModel::SampleParameters() {
	substitution_rate.sample();
}

void SingleSubstitutionRateModel::RecordState() {
	*substitution_model_out << substitution_rate << std::endl;
}

double SingleSubstitutionRateModel::SubstitutionProbability(int ancestral_state,
		int descendent_state, int site, double branch_length) {
	double probability = 0;

	if (ancestral_state != descendent_state)
		probability = 1 - std::exp(-substitution_rate * branch_length);
	else
		probability = std::exp(-substitution_rate * branch_length);

	return probability;
}

#include "SubstitutionRateMatrixModel.h"

#include <iostream>
#include <cmath> // For exp()
extern double Random();

int SubstitutionRateMatrixModel::number_of_substitution_rate_matrix_models = 0;

SubstitutionRateMatrixModel::SubstitutionRateMatrixModel() {
	id = number_of_substitution_rate_matrix_models;
	number_of_substitution_rate_matrix_models++;

	// Not necessary because the default constructor is called automatically
//	substitution_probability_matrix = Matrix<double>();
}

SubstitutionRateMatrixModel::SubstitutionRateMatrixModel(
		const SubstitutionRateMatrixModel& substitution_rate_matrix_model) {
	id = number_of_substitution_rate_matrix_models;
	number_of_substitution_rate_matrix_models++;

	substitution_rate_matrix = substitution_rate_matrix_model.substitution_rate_matrix;
	id = substitution_rate_matrix_model.id;
	substitution_model_out = substitution_rate_matrix_model.substitution_model_out;
}

SubstitutionRateMatrixModel::~SubstitutionRateMatrixModel() {
//	std::cout << "Single probability model deconstructor" << std::endl;
//	delete substitution_model_out;
}

SubstitutionRateMatrixModel* SubstitutionRateMatrixModel::Clone() {
//	std::cout << "cloning single probability model" << std::endl;
	return new SubstitutionRateMatrixModel(*this);
}

void SubstitutionRateMatrixModel::Initialize(
		int number_of_sites, vector<string> states) {
//	std::cout << "Initializing Single Probability Model" << std::endl;
	SubstitutionModel::Initialize(); // sets is_constant and is_reversible

	InitializeOutputStream(states);
	InitializeSubstitutionRateMatrix(states.size());
}

void SubstitutionRateMatrixModel::InitializeOutputStream(
		vector<string> states) {
	//	std::cout << "Initializing Single Probability Model" << std::endl;
	std::string substitution_model_out_file =
			"debug/substitution_rate_matrix_model_" + IdToString();
	substitution_model_out = new std::ofstream(
			substitution_model_out_file.c_str());

	// I do this for the probability matrix model. Perhaps it should be
	// abstracted away. Maybe it should inherit from a matrix substiution model
	// class that has the mechanics of printing a header from a vector of
	// states.
	int number_of_states = states.size();

	// I could print matrices two ways: either as a matrix or as a single line
	// I'm choosing a single line for now.
	// For now I am printing the diagonals too
	for (int row = 0; row < number_of_states; row++) {
		for (int column = 0; column < number_of_states; column++) {
			// if (row == column) continue; // To not print the diagonals
			if (not (row == 0 and column == 0))
				*substitution_model_out << "\t";

			*substitution_model_out << states.at(row)
					<< states.at(column);

		}
	}
	*substitution_model_out << std::endl;
}

void SubstitutionRateMatrixModel::InitializeSubstitutionRateMatrix(
		int number_of_states) {
	substitution_rate_matrix.resize(number_of_states, number_of_states, 0);
	SampleParameters();
}

/**
 * Notice that this is independent of the current substitution_probability.
 * This is the simplest (non-constant) way to sample a probability.
 * This samples from the uniform prior and not from the posterior or
 * marginal.
 *
 */
void SubstitutionRateMatrixModel::SampleParameters() {
	// These names could be row = state_from and column = state_to
	if (not is_constant) {
		for (int row = 0; row < substitution_rate_matrix.number_of_rows; row++) {
			double row_sum = 0;
			for (int column = 0;
					column < substitution_rate_matrix.number_of_columns; column++) {
				if (row == column)
					substitution_rate_matrix.at(row, column) = 0;
				else {
					// This is uniform prior, Metropolis-Hastings sampling
					// There also is a square uniform prior of 1 for values
					// from 0 to 1 and 0 after that
					// This sampling method forgets its previous state and
					// probably mixes poorly.
					substitution_rate_matrix.at(row, column) = Random();
					row_sum += substitution_rate_matrix.at(row, column);
				}
			}
			substitution_rate_matrix.at(row, row) = - row_sum;
		}
	}
}

void SubstitutionRateMatrixModel::RecordState() {
	for (int row = 0; row < substitution_rate_matrix.number_of_rows; row++) {
		for (int column = 0; column < substitution_rate_matrix.number_of_columns;
				column++) {
			// if (row == column) continue; // To not print the diagonals
			if (not (row == 0 and column == 0))
				*substitution_model_out << "\t";

			*substitution_model_out << substitution_rate_matrix.at(row, column);
		}
	}
	*substitution_model_out << std::endl;
}

double SubstitutionRateMatrixModel::SubstitutionProbability(int ancestral_state,
		int descendent_state, int site, double branch_length) {
	double probability = 0;

	if (ancestral_state != descendent_state)
		probability = 1
				- std::exp(
						-substitution_rate_matrix.at(ancestral_state,
								descendent_state) * branch_length);
	else
		probability = std::exp(
				-substitution_rate_matrix.at(ancestral_state, descendent_state)
						* branch_length);

	return probability;
}

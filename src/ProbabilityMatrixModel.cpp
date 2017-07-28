#include "ProbabilityMatrixModel.h"

#include <iostream>
#include <fstream>
#include <algorithm>

//#include algorithm

extern double Random();

int ProbabilityMatrixModel::number_of_probability_matrix_models = 0;

ProbabilityMatrixModel::ProbabilityMatrixModel() {
	id = number_of_probability_matrix_models;
	number_of_probability_matrix_models++;

	// Not necessary because the default constructor is called automatically
//	substitution_probability_matrix = Matrix<double>();
}

ProbabilityMatrixModel::ProbabilityMatrixModel(
		const ProbabilityMatrixModel& probability_matrix_model) {
	id = number_of_probability_matrix_models;
	number_of_probability_matrix_models++;

	substitution_probability_matrix =
			probability_matrix_model.substitution_probability_matrix;
	id = probability_matrix_model.id;
	substitution_model_out = probability_matrix_model.substitution_model_out;
}

ProbabilityMatrixModel::~ProbabilityMatrixModel() {
//	std::cout << "Single probability model deconstructor" << std::endl;
//	delete substitution_model_out;
}

ProbabilityMatrixModel* ProbabilityMatrixModel::Clone() {
//	std::cout << "cloning single probability model" << std::endl;
	return new ProbabilityMatrixModel(*this);
}

void ProbabilityMatrixModel::Initialize(int number_of_sites,
		vector<string> states) {
//	std::cout << "Initializing Single Probability Model" << std::endl;
	SubstitutionModel::Initialize(); // sets is_constant


	InitializeStateFromFile("matrix", states);
	// Will call InitializeState(number_of_sites,
	InitializeProbabilityMatrix(states.size());
	InitializeOutputStream(states);
}

void ProbabilityMatrixModel::InitializeOutputStream(vector<string> states) {
	//	std::cout << "Initializing Single Probability Model" << std::endl;
	std::string substitution_model_out_file = "debug/probability_matrix_model_"
			+ IdToString();
	substitution_model_out = new std::ofstream(
			substitution_model_out_file.c_str());

	int number_of_states = states.size();

	// I could print matrices two ways: either as a matrix or as a single line
	// I'm choosing a single line for now.
	// For now I am printing the diagonals too
	for (int row = 0; row < number_of_states; row++) {
		for (int column = 0; column < number_of_states; column++) {
			// if (row == column) continue; // To not print the diagonals
			if (not (row == 0 and column == 0))
				*substitution_model_out << "\t";

			*substitution_model_out << states.at(row) << states.at(column);

		}
	}
	*substitution_model_out << std::endl;
}

void ProbabilityMatrixModel::InitializeProbabilityMatrix(int number_of_states) {
	substitution_probability_matrix.resize(number_of_states, number_of_states,
			0);
	SampleParameters();
}

/**
 * Notice that this is independent of the current substitution_probability.
 * This is the simplest (non-constant) way to sample a probability.
 * This samples from the uniform prior and not from the posterior or
 * marginal.
 *
 */
void ProbabilityMatrixModel::SampleParameters() {
	// These names could be row = state_from and column = state_to
	if (not is_constant) {
		for (int row = 0; row < substitution_probability_matrix.number_of_rows;
				row++) {
			double row_sum = 0;
			for (int column = 0;
					column < substitution_probability_matrix.number_of_columns;
					column++) {
				if (row == column)
					substitution_probability_matrix.at(row, column) = 0;
				else {
					// This is uniform prior, Metropolis-Hastings sampling
					substitution_probability_matrix.at(row, column) = Random();
					row_sum += substitution_probability_matrix.at(row, column);
				}
			}

			for (int column = 0;
					column < substitution_probability_matrix.number_of_columns;
					column++) {
				substitution_probability_matrix.at(row, column) /= row_sum;
			}
		}
	}
}

void ProbabilityMatrixModel::RecordState() {
	for (int row = 0; row < substitution_probability_matrix.number_of_rows;
			row++) {
		for (int column = 0;
				column < substitution_probability_matrix.number_of_columns;
				column++) {
			// if (row == column) continue; // To not print the diagonals
			if (not (row == 0 and column == 0))
				*substitution_model_out << "\t";

			*substitution_model_out
					<< substitution_probability_matrix.at(row, column);
		}
	}
	*substitution_model_out << std::endl;
}

double ProbabilityMatrixModel::SubstitutionProbability(int ancestral_state,
		int descendent_state, int site, double branch_length) {
	double probability = 0;

	if (ancestral_state != descendent_state)
		probability = substitution_probability_matrix.at(ancestral_state,
				descendent_state);
	else
		probability = 1
				- substitution_probability_matrix.at(ancestral_state,
						descendent_state);

	return probability;
}

void ProbabilityMatrixModel::InitializeStateFromFile(
		std::string state_in_file, vector<string> states) {
	std::ifstream state_in(state_in_file.c_str());
	int number_of_pairs = states.size() * states.size();

	// read header
	vector<string> labels(number_of_pairs);
	for (int pair = 0; pair < number_of_pairs; pair++) {
		state_in >> labels.at(pair);
		std::cout << labels.at(pair) << " ";
	}
	std::cout << std::endl;

	Matrix<double> matrix(states.size(), states.size());
	for (int label = 0; label < number_of_pairs; label++) {
		// This is instead of index to state map
		int row = find(states.begin(), states.end(),
				labels.at(label).substr(0, 1)) - states.begin();
		int column = find(states.begin(), states.end(),
				labels.at(label).substr(1, 2)) - states.begin();
		std::cout << row << " " << column << std::endl;
		state_in >> substitution_probability_matrix.at(row, column);
	}

	for (int row = 0; row < matrix.number_of_rows; row++) {
		for (int column = 0; column < matrix.number_of_columns; column++) {
			std::cout << substitution_probability_matrix.at(row, column) << " ";
		}
		std::cout << std::endl;
	}
	exit(1);
}

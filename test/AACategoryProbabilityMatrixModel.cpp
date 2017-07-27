#include "AACategoryProbabilityMatrixModel.h"

#include <iostream>

extern double Random();

int AACategoryProbabilityMatrixModel::number_of_probability_matrix_models = 0;

AACategoryProbabilityMatrixModel::AACategoryProbabilityMatrixModel() {
	id = number_of_probability_matrix_models;
	number_of_probability_matrix_models++;

	// Not necessary because the default constructor is called automatically
//	substitution_probability_matrix = Matrix<double>();
}

AACategoryProbabilityMatrixModel::AACategoryProbabilityMatrixModel(
		const AACategoryProbabilityMatrixModel& probability_matrix_model) {
	id = number_of_probability_matrix_models;
	number_of_probability_matrix_models++;

	substitution_probability_matrix =
			probability_matrix_model.substitution_probability_matrix;
	id = probability_matrix_model.id;
	substitution_model_out = probability_matrix_model.substitution_model_out;
}

AACategoryProbabilityMatrixModel::~AACategoryProbabilityMatrixModel() {
//	std::cout << "Single probability model deconstructor" << std::endl;
//	delete substitution_model_out;
}

AACategoryProbabilityMatrixModel* AACategoryProbabilityMatrixModel::Clone() {
//	std::cout << "cloning single probability model" << std::endl;
	return new AACategoryProbabilityMatrixModel(*this);
}

void AACategoryProbabilityMatrixModel::Initialize(int number_of_sites,
		std::map<int, std::string> integer_to_state) {
//	std::cout << "Initializing Single Probability Model" << std::endl;
	SubstitutionModel::Initialize();

	InitializeCategories();

	int number_of_categories = 3; //Arbitrarily chosen. Instead read from file
	InitializeProbabilityMatrix(number_of_categories);

	InitializeOutputStream();
}

void AACategoryProbabilityMatrixModel::InitializeOutputStream(
		) {
	//	std::cout << "Initializing Single Probability Model" << std::endl;
	std::string substitution_model_out_file = "debug/AACategory_probability_matrix_model_"
			+ IdToString();
	substitution_model_out = new std::ofstream(
			substitution_model_out_file.c_str());

	int number_of_categories = residue_to_category.size();

	// I could print matrices two ways: either as a matrix or as a single line
	// I'm choosing a single line for now.
	// For now I am printing the diagonals too
	for (int row = 0; row < number_of_categories; row++) {
		for (int column = 0; column < number_of_categories; column++) {
			// if (row == column) continue; // To not print the diagonals
			if (not (row == 0 and column == 0))
				*substitution_model_out << "\t";

			*substitution_model_out << residue_to_category[row] << "_"
					<< residue_to_category[column];
		}
	}
	*substitution_model_out << std::endl;
}

void AACategoryProbabilityMatrixModel::InitializeCategories() {

	// Read in categories here instead of make them up...
	for (int i = 0; i < 20; i++) {
		if (i < 5) {
			residue_to_category[i] = 0;
		} else if (i < 14) {
			residue_to_category[i] = 1;
		} else {
			residue_to_category[i] = 2;
		}
	}
}
void AACategoryProbabilityMatrixModel::InitializeProbabilityMatrix(
		int number_of_categories) {
	substitution_probability_matrix.resize(number_of_categories,
			number_of_categories, 0);
	SampleParameters();
}

/**
 * Notice that this is independent of the current substitution_probability.
 * This is the simplest (non-constant) way to sample a probability.
 * This samples from the uniform prior and not from the posterior or
 * marginal.
 *
 */
void AACategoryProbabilityMatrixModel::SampleParameters() {
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

void AACategoryProbabilityMatrixModel::RecordState() {
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

double AACategoryProbabilityMatrixModel::SubstitutionProbability(
		int ancestral_state, int descendent_state, int site,
		double branch_length) {
	double probability = 0;

	// Notice how the model transforms the residue into the category
	if (ancestral_state != descendent_state)
		probability = substitution_probability_matrix.at(
				residue_to_category[ancestral_state],
				residue_to_category[descendent_state]);
	else
		probability = 1
				- substitution_probability_matrix.at(
						residue_to_category[ancestral_state],
						residue_to_category[descendent_state]);

	return probability;
}

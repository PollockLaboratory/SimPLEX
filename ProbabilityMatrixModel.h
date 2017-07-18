#ifndef ProbabilityMatrixModel_h_
#define ProbabilityMatrixModel_h_

#include <fstream>
#include <map>

#include "SubstitutionModel.h"

#include "Matrix.h"

/**
 * This substitution model has a probability for a substitution between an
 * ancestral state and a descendant state.
 *
 * T
 *
 */

class ProbabilityMatrixModel: public SubstitutionModel {

public:
	ProbabilityMatrixModel();
	ProbabilityMatrixModel(const ProbabilityMatrixModel&);
	virtual ~ProbabilityMatrixModel();

	virtual ProbabilityMatrixModel* Clone();

	virtual void Initialize(int number_of_state, vector<string> states);

	virtual void SampleParameters();

	virtual void RecordState();

	virtual double SubstitutionProbability(int ancestral_state,
			int descendent_state, int site, double branch_length);
private:
	static int number_of_probability_matrix_models;
	Matrix<double> substitution_probability_matrix;

	virtual void InitializeOutputStream(vector<string> states);
	void InitializeProbabilityMatrix(int number_of_states);

	virtual void InitializeStateFromFile(std::string state_in_file,
			vector<string> states);

	virtual void InitializeStateFromFile(std::string state_in_file) {}

};

#endif

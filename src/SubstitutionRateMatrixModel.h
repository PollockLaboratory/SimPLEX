#ifndef SubstitutionRateMatrixModel_h_
#define SubstitutionRateMatrixModel_h_

#include <fstream>
#include <map>

#include "SubstitutionModel.h"

#include "Matrix.h"

/**
 * This substitution model has a substitution rate for a substitution between an
 * ancestral state and a descendant state.
 *
 * Stephen Pollard
 * 12/22/2013
 *
 *
 */

class SubstitutionRateMatrixModel: public SubstitutionModel {

public:
	SubstitutionRateMatrixModel();
	SubstitutionRateMatrixModel(const SubstitutionRateMatrixModel&);
	virtual ~SubstitutionRateMatrixModel();

	virtual SubstitutionRateMatrixModel* Clone();

	virtual void Initialize(int number_of_states, vector<string> states);

	virtual void SampleParameters();

	virtual void RecordState();

	virtual double SubstitutionProbability(int ancestral_state,
			int descendent_state, int site, double branch_length);
private:
	static int number_of_substitution_rate_matrix_models;
	Matrix<double> substitution_rate_matrix;

	virtual void InitializeOutputStream(vector<string> states);
	void InitializeSubstitutionRateMatrix(int number_of_states);

	virtual void InitializeStateFromFile(std::string) {};


};

#endif

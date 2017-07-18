#ifndef AACategoryProbabilityMatrixModel_h_
#define AACategoryProbabilityMatrixModel_h_

#include <fstream>
#include <map>

#include "SubstitutionModel.h"

#include "Matrix.h"

/**
 * Still need to read in the categories from a file.
 *
 */

class AACategoryProbabilityMatrixModel: public SubstitutionModel {

public:
	AACategoryProbabilityMatrixModel();
	AACategoryProbabilityMatrixModel(const AACategoryProbabilityMatrixModel&);
	virtual ~AACategoryProbabilityMatrixModel();

	virtual AACategoryProbabilityMatrixModel* Clone();

	virtual void Initialize(int number_of_sites,
			std::map<int, std::string> integer_to_state);

	virtual void SampleParameters();

	virtual void RecordState();

	virtual double SubstitutionProbability(int ancestral_state,
			int descendent_state, int site, double branch_length);
private:
	static int number_of_probability_matrix_models;
	Matrix<double> substitution_probability_matrix;
	std::map<int, int> residue_to_category;


	void InitializeCategories();
	void InitializeProbabilityMatrix(int number_of_states);
	void InitializeOutputStream();
};

#endif

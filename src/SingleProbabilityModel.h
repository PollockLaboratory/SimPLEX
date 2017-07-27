#ifndef SingleProbabilityModel_h_
#define SingleProbabilityModel_h_

#include <fstream>

#include "SubstitutionModel.h"

/**
 * This substitution model has a single probability for a substitution between
 * any two residues along a branch of any length. The probability of not seeing
 * a substitution is 1 minus the probability of seeing a substitution.
 *
 * The substitution probability is not constant; I think the best way to make
 * a model constant is to never call SampleParameters() on it.
 *
 */

class SingleProbabilityModel: public SubstitutionModel {

public:
	SingleProbabilityModel();
	SingleProbabilityModel(const SingleProbabilityModel&);
	virtual ~SingleProbabilityModel();

	virtual SingleProbabilityModel* Clone();

	virtual void Initialize(int number_of_sites, std::vector<std::string> states);

	virtual void SampleParameters();

	virtual void RecordState();

	virtual double SubstitutionProbability(int ancestral_state,
			int descendent_state, int site, double branch_length);
private:
	static int number_of_single_probability_models;
	// Might want to call this substitution_probability like substitution_rate
	// or vice versa.
	double substitution_probability;

	void InitializeOutputStream();
	virtual void InitializeStateFromFile(std::string state_in_file);
};

#endif

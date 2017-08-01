#ifndef SingleSubstitutionRateModel_h_
#define SingleSubstitutionRateModel_h_

#include <fstream>

#include "SubstitutionModel.h"

#include "FloatParameter.h"

/**
 * Single substitution rate model
 *
 * A better name for it is a Substitution rate model because it is specifically
 * dealing with substitution substitutions and not substitutions in general. For
 * example, it does not account for insertions or deletions which are substitutions
 * but not substitutions. Therefore a better name is substitution rate model.
 *
 */

class SingleSubstitutionRateModel: public SubstitutionModel {

public:
	SingleSubstitutionRateModel();
	virtual ~SingleSubstitutionRateModel();

	virtual SingleSubstitutionRateModel* Clone();

	virtual void Initialize(int number_of_sites, std::vector<std::string> states);
	virtual void InitializeStateFromFile(std::string state_in_file);

	virtual void SampleParameters();

	virtual void RecordState();

	virtual double SubstitutionProbability(int ancestral_state,
			int descendent_state, int site, double branch_length);
private:
	static int number_of_single_substitution_rate_models;
	FloatParameter substitution_rate;

	void InitializeOutputStream();
};

#endif

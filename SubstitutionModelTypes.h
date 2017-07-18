#ifndef SubstitutionModelTypes_h_
#define SubstitutionModelTypes_h_

#include "SubstitutionModel.h"
#include "SingleProbabilityModel.h"
#include "SingleSubstitutionRateModel.h"
#include "ProbabilityMatrixModel.h"
#include "SubstitutionRateMatrixModel.h"

#include "RadiusDependentModel.h"

#include "SubstitutionMixtureModel.h"

#include "Options.h"
extern Options options;

/**
 * I'm not sure of the best way to do this. These are not globals; they only
 * are needed in two files: Model.cpp and SubsitutionMixtureModel.cpp .
 *
 * And now RadiusDependentModel.cpp.
 *
 * I might use enum's?
 *
 */

static const int single_probability_model_type = 0;
static const int single_substitution_rate_model_type = 1;
static const int probability_matrix_model_type = 2;
static const int substitution_rate_matrix_model_type = 3;
static const int radius_dependent_model_type = 4;
static const int substitution_models_mixture_model_type = 5;

static SubstitutionModel* GetNextSubstitutionModel() {
	SubstitutionModel* substitution_model = NULL;

	// No, I will not allow recycling of model types. They must be specified
		// correctly or fail.
	if (options.substitution_model_types.size() == 0) {
		std::cerr << "Not enough substitution model types given" << std::endl;
		exit(-1);
	}

	int next_model_type = options.substitution_model_types.front();
	options.substitution_model_types.pop();

	// This is a good candidate for a switch or case construction
	if (next_model_type == single_probability_model_type) {
		substitution_model = new SingleProbabilityModel();
	} else if (next_model_type == single_substitution_rate_model_type) {
		substitution_model = new SingleSubstitutionRateModel();
	} else if (next_model_type == probability_matrix_model_type) {
		substitution_model = new ProbabilityMatrixModel();
	} else if (next_model_type == substitution_rate_matrix_model_type) {
		substitution_model = new SubstitutionRateMatrixModel();
	} else if (next_model_type == radius_dependent_model_type) {
		substitution_model = new RadiusDependentModel();
	} else if (next_model_type == substitution_models_mixture_model_type) {
		substitution_model = new SubstitutionMixtureModel();
	} else {
		std::cerr
				<< "Substitution model type not recognized in get next sub model "
				<< next_model_type << std::endl;
		std::exit(-1);
	}
	return substitution_model;
}

#endif

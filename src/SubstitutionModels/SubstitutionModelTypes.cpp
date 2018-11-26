#include <iostream>
#include <cstdlib>

#include "SubstitutionModel.h"
#include "../Environment.h"

#include "Types/GeneralTimeReversible.h"
#include "Types/SingleRate.h"

extern Environment env;

static const int general_time_reversible_model_type = 0;
static const int single_rate_model_type = 1;

SubstitutionModel* GetSubstitutionModel() {
	/*
	 * Given the global Options will dynamically attach the desired substitution model to the 
	 * substitution_model pointer. The model remains un initialized.
	 */

	SubstitutionModel* substitution_model = NULL;
	int next_model_type = env.get_int("substitution_model_type");

	if (next_model_type == general_time_reversible_model_type) {
		substitution_model = new GeneralTimeReversible();
	} else if (next_model_type == single_rate_model_type) {
		substitution_model = new SingleRate();
	} else {
		std::cerr << "Substitution model type not recognized in get next sub model " << next_model_type << std::endl;
		std::exit(-1);
	}

	return substitution_model;
}

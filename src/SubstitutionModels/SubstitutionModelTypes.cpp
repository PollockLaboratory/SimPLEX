#include <iostream>
#include <cstdlib>

#include "SubstitutionModel.h"
#include "../Environment.h"

#include "Types/GeneralTimeReversible.h"
#include "Types/CustomModel.h"

extern Environment env;

// This should be ENUM;
static const int general_time_reversible_model_type = 0;
// static const int single_rate_model_type = 1;
// static const int fixed_rate_model_type = 2;
// static const int category_GTR_model_type = 3;
static const int custom_model = 1;

SubstitutionModel* GetSubstitutionModel() {
  /*
   * Given the global Options will dynamically attach the desired substitution model to the 
   * substitution_model pointer. The model remains un initialized.
   */

  SubstitutionModel* substitution_model = NULL;
  int next_model_type = env.get_int("substitution_model_type");

  if (next_model_type == general_time_reversible_model_type) {
    substitution_model = new GeneralTimeReversible();
  } else if (next_model_type == custom_model) {
    substitution_model = new CustomModel();
  } else {
    std::cerr << "Substitution model type not recognized in get next sub model " << next_model_type << std::endl;
    std::exit(EXIT_FAILURE);
  }

  return substitution_model;
}

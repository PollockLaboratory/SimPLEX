#ifndef GeneralTimeReversible_h_
#define GeneralTimeReversible_h_

#include <fstream>

#include "SubstitutionModel.h"
#include "../../Parameters/ParameterSet.h"
#include "../../Parameters/Types/ContinuousFloat.h"

class GeneralTimeReversible: public SubstitutionModel {
	public:
		GeneralTimeReversible();
		virtual void Initialize(int number_of_sites, std::vector<std::string> states);
		virtual double SubstitutionProbability(int ancestral_state, int descendent_state, int site, double branch_length);
	private:
};

#endif

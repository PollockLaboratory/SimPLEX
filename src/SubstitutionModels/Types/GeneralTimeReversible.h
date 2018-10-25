#ifndef GeneralTimeReversible_h_
#define GeneralTimeReversible_h_

#include <fstream>

#include "SubstitutionModel.h"
#include "../../Parameters/ParameterSet.h"
#include "../../Parameters/Types/ContinuousFloat.h"
#include "../../Parameters/Types/VirtualSubstitutionRate.h"
#include "../../Parameters/AbstractValue.h"
#include "../../Parameters/RateVector.h"

class GeneralTimeReversible: public SubstitutionModel {
	public:
		GeneralTimeReversible();
		virtual void Initialize(int number_of_sites, std::vector<std::string> states);
	private:
};

#endif

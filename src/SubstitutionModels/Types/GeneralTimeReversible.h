#ifndef GeneralTimeReversible_h_
#define GeneralTimeReversible_h_

#include <fstream>

#include "SubstitutionModel.h"
#include "../Components/ComponentSet.h"

#include "../Components/Types/AbstractValueTypes.h"
#include "../Components/Types/SampleableValueTypes.h"

#include "../Components/AbstractComponent.h"
#include "../Components/RateVector.h"

class GeneralTimeReversible: public SubstitutionModel {
	public:
		GeneralTimeReversible();
		virtual void Initialize(int number_of_sites, std::vector<std::string> states);
	private:
};

#endif

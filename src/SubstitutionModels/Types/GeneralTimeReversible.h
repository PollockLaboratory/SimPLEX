#ifndef GeneralTimeReversible_h_
#define GeneralTimeReversible_h_

#include <fstream>

#include "SubstitutionModel.h"
#include "../Components/ComponentSet.h"

#include "../Components/AbstractValueTypes.h"
#include "../Components/SampleableValueTypes.h"

#include "../Components/AbstractComponent.h"
#include "../Components/RateVector.h"

class GeneralTimeReversible: public SubstitutionModel {
	public:
		GeneralTimeReversible();
		virtual void Initialize();
		virtual void Initialize(int number_of_sites, std::vector<std::string> states);
	private:
};

#endif

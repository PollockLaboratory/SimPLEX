#ifndef SingleRate_h_
#define SingleRate_h_

#include <fstream>

#include "SubstitutionModel.h"
#include "../Components/ComponentSet.h"

#include "../Components/Types/AbstractValueTypes.h"
#include "../Components/Types/SampleableValueTypes.h"

#include "../Components/AbstractComponent.h"
#include "../Components/RateVector.h"

class SingleRate: public SubstitutionModel {
	public:
		SingleRate();
		virtual void Initialize(int number_of_sites, std::vector<std::string> states);
	private:
};

#endif

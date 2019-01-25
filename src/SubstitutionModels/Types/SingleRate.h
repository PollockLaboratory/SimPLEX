#ifndef SingleRate_h_
#define SingleRate_h_

#include <fstream>

#include "SubstitutionModel.h"
#include "../../Parameters/ComponentSet.h"
#include "../../Parameters/Types/ContinuousFloat.h"
#include "../../Parameters/Types/VirtualSubstitutionRate.h"
#include "../../Parameters/AbstractValue.h"
#include "../../Parameters/RateVector.h"

class SingleRate: public SubstitutionModel {
	public:
		SingleRate();
		virtual void Initialize(int number_of_sites, std::vector<std::string> states);
	private:
};

#endif

#ifndef CategoryGTR_h
#define CategoryGTR_h

#include <fstream>

#include "SubstitutionModel.h"
#include "../Components/ComponentSet.h"

#include "../Components/Types/AbstractValueTypes.h"
#include "../Components/Types/SampleableValueTypes.h"
#include "../Components/ComponentTypes.h"

#include "../Components/AbstractComponent.h"
#include "../Components/RateVector.h"

class CategoryGTR: public SubstitutionModel {
	public:
		CategoryGTR();
		virtual void Initialize(int number_of_sites, std::vector<std::string> states);
	private:
};

#endif

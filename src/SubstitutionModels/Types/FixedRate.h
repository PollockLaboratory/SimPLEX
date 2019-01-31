#ifndef FixedRate_h_
#define FixedRate_h_

#include <fstream>

#include "SubstitutionModel.h"
#include "../Components/ComponentSet.h"
#include "../Components/Types/AbstractValueTypes.h"
#include "../Components/AbstractComponent.h"
#include "../Components/RateVector.h"

class FixedRate: public SubstitutionModel {
 public:
  FixedRate();
  virtual void Initialize(int number_of_sites, std::vector<std::string> states);
 private:
};

#endif

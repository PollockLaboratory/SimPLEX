#ifndef FixedRate_h_
#define FixedRate_h_

#include <fstream>

#include "SubstitutionModel.h"
#include "../../Parameters/ParameterSet.h"
#include "../../Parameters/Types/FixedFloat.h"
#include "../../Parameters/Types/VirtualSubstitutionRate.h"
#include "../../Parameters/AbstractValue.h"
#include "../../Parameters/RateVector.h"

class FixedRate: public SubstitutionModel {
 public:
  FixedRate();
  virtual void Initialize(int number_of_sites, std::vector<std::string> states);
 private:
};

#endif

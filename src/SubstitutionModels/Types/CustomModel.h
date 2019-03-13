#ifndef CustomModel_h
#define CustomModel_h

#include <fstream>
#include <map>

#include "SubstitutionModel.h"
#include "../Components/ComponentSet.h"

#include "../Components/AbstractValueTypes.h"
#include "../Components/SampleableValueTypes.h"
#include "../Components/ComponentTypes.h"

#include "../Components/AbstractComponent.h"
#include "../Components/RateVector.h"

#include "sol2/sol.hpp"

class CustomModel: public SubstitutionModel {
 public:
  CustomModel();

  virtual void Initialize();
  virtual void Initialize(int number_of_sites, std::vector<std::string> states);

  void set_name(std::string);
  void set_states(sol::table);
 private:
  std::string name;
};

#endif

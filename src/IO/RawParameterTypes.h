#ifndef RawParameterTypes_h_
#define RawParameterTypes_h_

#include <iostream>
#include <list>

#include "../ModelParts/SubstitutionModels/Parameters.h"
#include "sol2/sol.hpp"

namespace IO { 
  class ParameterWrapper {
  private:
    std::string name;
  public:
    ParameterWrapper(AbstractComponent*);
    std::string get_name();
    std::string get_type();
    AbstractComponent* parameter;
  };

  ParameterWrapper new_parameter(std::string, std::string, sol::table);
  ParameterWrapper new_categories(std::string, sol::table);
}

#endif

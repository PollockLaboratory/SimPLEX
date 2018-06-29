#include <iostream>

#include "AbstractParameter.h"

AbstractParameter::AbstractParameter(std::string parameter_name) {
	name = parameter_name;
	fixedQ = true;
}



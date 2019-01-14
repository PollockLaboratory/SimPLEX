#include <iostream>
#include "VirtualSubstitutionRate.h"

VirtualSubstitutionRate::VirtualSubstitutionRate(std::string parameter_name, double unif) : AbstractValue(parameter_name) {
	u = unif;
}

double VirtualSubstitutionRate::getValue() {
	return(value);
}

double  VirtualSubstitutionRate::getOldValue() {
	return(previous_value);
}

void VirtualSubstitutionRate::printValue() {
	std::cout << "VirtualSubstitutionRate - " << name << ": " << value << std::endl;
}

void VirtualSubstitutionRate::refresh() {
  previous_value = value;
  double total = 0.0;
  for(auto it = dependent_values.begin(); it != dependent_values.end(); ++it) {
    total += (*it)->getValue();
  }

  value = u - total;
  if(value < 0.0 || value > 1.0) {
    throw OutOfBoundsException("VirtualSubstitutionRate out of bounds.");
  }
}


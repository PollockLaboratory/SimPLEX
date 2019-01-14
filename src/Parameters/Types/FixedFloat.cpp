#include <iostream>
#include "FixedFloat.h"

FixedFloat::FixedFloat(std::string parameter_name, double v) : AbstractValue(parameter_name) {
  value = v;
}

double FixedFloat::getValue() {
  return(value);
}

double  FixedFloat::getOldValue() {
  return(value);
}

void FixedFloat::printValue() {
  std::cout << "FixedFloat - " << name << ": " << value << std::endl;
}

void FixedFloat::refresh() {	
}

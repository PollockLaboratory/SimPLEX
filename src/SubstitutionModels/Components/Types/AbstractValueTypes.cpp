#include "AbstractValueTypes.h"

// FIXED FLOAT

FixedFloat::FixedFloat(std::string parameter_name, double v) : AbstractValue(parameter_name) {
  value = v;
}

double FixedFloat::getValue() {
  return(value);
}

double  FixedFloat::getOldValue() {
  return(value);
}

void FixedFloat::print() {
  std::cout << "FixedFloat - " << name << ": " << value << std::endl;
}

void FixedFloat::refresh() {	
}

// VIRTUAL SUBSTITUTION RATE

VirtualSubstitutionRate::VirtualSubstitutionRate(std::string parameter_name, double unif) : AbstractValue(parameter_name) {
	u = unif;
}

double VirtualSubstitutionRate::getValue() {
	return(value);
}

double  VirtualSubstitutionRate::getOldValue() {
	return(previous_value);
}

void VirtualSubstitutionRate::print() {
	std::cout << "VirtualSubstitutionRate - " << name << ": " << value << std::endl;
}

void VirtualSubstitutionRate::refresh() {
  previous_value = value;
  double total = 0.0;
  for(auto it = dependent_rates.begin(); it != dependent_rates.end(); ++it) {
    total += (*it)->getValue();
  }

  value = u - total;
  if(value < 0.0 || value > 1.0) {
    throw OutOfBoundsException("VirtualSubstitutionRate out of bounds.");
  }
}

void VirtualSubstitutionRate::add_rate(AbstractValue* v) {
  add_dependancy(v);
  dependent_rates.push_back(v);
}



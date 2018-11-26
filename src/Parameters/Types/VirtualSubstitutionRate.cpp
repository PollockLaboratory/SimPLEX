#include <iostream>
#include "VirtualSubstitutionRate.h"

VirtualSubstitutionRate::VirtualSubstitutionRate(std::string parameter_name, double unif) : AbstractDependentParameter(parameter_name) {
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
}

void VirtualSubstitutionRate::add_dependancy(AbstractValue* v) {
	dependent_values.push_back(v);
}

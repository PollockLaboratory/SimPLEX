#include <iostream>
#include <stdlib.h> //This gives rand.

#include "ContinuousFloat.h" 

ContinuousFloat::ContinuousFloat(std::string name, double initial_value = 0.0, double initial_std_dev = 1.0) : AbstractParameter(name) {
	/*
	 * The default constructor for the Continuous Float parameter class.
	 */
	value = initial_value;
	std_dev = initial_std_dev;
	previous_value = 0.0;
}

void ContinuousFloat::printValue() {
	std::cout << "Continuous float - " << name << ": " << value << std::endl;
}

void ContinuousFloat::sample() {
	previous_value = value; 
	fixedQ = false;

	double r = ((rand() % 10000) / 10000.0) - 0.5;
	value = value + (r * std_dev);
}

double ContinuousFloat::getValue() {
	return value;
}

void ContinuousFloat::undo() {
	value = previous_value;
	previous_value = 0;
	fixedQ = true;
}

void ContinuousFloat::fix() {
	previous_value = 0;
	fixedQ = true;
}


#include <iostream>
#include <stdlib.h> //This gives rand.
#include <limits>

#include "ContinuousFloat.h" 

double inf = std::numeric_limits<double>::infinity();

ContinuousFloat::ContinuousFloat(std::string name, double initial_value = 0.0, double initial_std_dev = 1.0) : AbstractParameter(name), value(initial_value), std_dev(initial_std_dev) {
	/*
	 * The default constructor for the Continuous Float parameter class.
	 */
	lower_bound = -inf;
	upper_bound = inf;

	previous_value = 0.0;
}

ContinuousFloat::ContinuousFloat(std::string name, double initial_value = 0.0, double initial_std_dev = 1.0, double lower_bound = -inf) : AbstractParameter(name), value(initial_value), std_dev(initial_std_dev), lower_bound(lower_bound) {
	/*
	 * The default constructor for the Continuous Float parameter class - with bounds
	 */
	upper_bound = inf;

	previous_value = 0.0;
}

ContinuousFloat::ContinuousFloat(std::string name, double initial_value = 0.0, double initial_std_dev = 1.0, double lower_bound = -inf, double upper_bound = inf) : AbstractParameter(name), value(initial_value), std_dev(initial_std_dev), lower_bound(lower_bound), upper_bound(upper_bound) {
	/*
	 * The default constructor for the Continuous Float parameter class - with bounds
	 */
	previous_value = 0.0;
}

// Utils
void ContinuousFloat::printValue() {
	std::cout << "Continuous float - " << name << ": " << value << std::endl;
}

bool ContinuousFloat::sample() {
	previous_value = value; 
	fixedQ = false;

	double r = ((rand() % 10000) / 10000.0) - 0.5;
	value = value + (r * std_dev);

	if(value < lower_bound) {
		value = 2*lower_bound - value;
	}

	if(value > upper_bound) {
		value = 2*upper_bound - value;
	}
	return(true);
}

double ContinuousFloat::getValue() {
	return(value);
}

double ContinuousFloat::getOldValue() {
	if(fixedQ) {
		std::cout << "Error: in ContinuousFloat::getOldValue - already fixed." << std::endl;
		exit(EXIT_FAILURE);
	}
	return(previous_value);
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


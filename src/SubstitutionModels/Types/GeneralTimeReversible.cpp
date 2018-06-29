#include "GeneralTimeReversible.h"

#include <iostream>
#include <cstdlib>

#include "Options.h"

extern Options options;

GeneralTimeReversible::GeneralTimeReversible() {
}

void GeneralTimeReversible::Initialize(int number_of_sites, std::vector<std::string> states) {
	/*
	 * std::cout << "Initializing Single Probability Model" << std::endl;
	 */
	//Hamish	
	substitution_model_out = CreateOutputStream("Single_probability_model");

	std::cout << "Initializing General Time Reversible Model (GTR)." << std::endl;

	parameters.add_parameter(new ContinuousFloat("a", 0.5, 0.1));
	parameters.add_parameter(new ContinuousFloat("b", 0.5, 0.1));
	parameters.add_parameter(new ContinuousFloat("c", 0.6, 0.1));
	parameters.print();
	parameters.Initialize(substitution_model_out);
}

double GeneralTimeReversible::SubstitutionProbability(int ancestral_state, int descendent_state, int site, double branch_length) {
	/*
	 * Gives the likelihood of a particular sustitution.
	 */
	//std::cout << "Calculating probability " << ancestral_state << " " << descendent_state << " pos " << site << " branch length" << std::endl;

	double probability = 0;
	double rate; 

	if(ancestral_state < 10) { 
		rate = parameters.get("a");
	} else {
		rate = parameters.get("b");
	}
	
	if(ancestral_state != descendent_state) {
		probability = rate;
	} else {
		probability = 1.0 - rate;
	}

	return probability;
}

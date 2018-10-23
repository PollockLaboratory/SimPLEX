#include "GeneralTimeReversible.h"

#include <iostream>
#include <cstdlib>

#include "Environment.h"

extern Environment env;

GeneralTimeReversible::GeneralTimeReversible() {
}

void GeneralTimeReversible::Initialize(int number_of_sites, std::vector<std::string> states) {
	/*
	 * std::cout << "Initializing Single Probability Model" << std::endl;
	 */

	std::cout << "Initializing General Time Reversible Model (GTR)." << std::endl;

	std::array<std::string, 20> aa = {"A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"};
	
	add_rate_matrix(new RateMatrix("Q", 20));
	finalize();
}


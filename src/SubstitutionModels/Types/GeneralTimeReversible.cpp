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

	float u = env.u;

	std::cout << "Protein General Time Reversible Model (GTR)." << std::endl;

	std::array<std::string, 20> aa = {"A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"};

	std::vector<std::vector<AbstractValue*>> Q(20, std::vector<AbstractValue*>(20, NULL));

	AbstractValue* r = NULL;
	for(int i = 0; i < 20; i++) {
		for(int j = 0; j < 20; j++) {
			if(i == j) {
				Q[i][j] = new VirtualSubstitutionRate(aa[i] + aa[j], u);
			} else if(i < j) {
				Q[i][j] = new ContinuousFloat(aa[i] + aa[j], 0.01, 0.0001, 0.0);
			} else {
				Q[i][j] = Q[j][i];
			}
		}

		for(int j = 0; j < 20; j++) {
			if(i != j) {
				dynamic_cast<AbstractDependentParameter*>(Q[i][i])->add_dependancy(Q[i][j]);
			}
		}

		RateVector* v = new RateVector("rv-" + aa[i], i, Q[i]);
		add_rate_vector(v);
	}
	
	finalize();
}

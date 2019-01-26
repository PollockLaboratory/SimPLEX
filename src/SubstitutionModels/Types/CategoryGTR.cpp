#include "CategoryGTR.h"

#include <iostream>
#include <cstdlib>

#include "Environment.h"

extern Environment env;

CategoryGTR::CategoryGTR() {
}

void CategoryGTR::Initialize(int number_of_sites, std::vector<std::string> states) {
	/*
	 * std::cout << "Initializing Single Probability Model" << std::endl;
	 */

  float u = env.u;

  std::cout << "Protein General Time Reversible Model (GTR)." << std::endl;

  std::array<std::string, 20> aa = {"A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"};

  std::vector<std::vector<AbstractValue*>> Q(20, std::vector<AbstractValue*>(20, NULL));

  RateCategories* rc = new RateCategories("Rate Categories", 0.0, 0.1, 100);

  AbstractValue* r = NULL;
  for(int i = 0; i < 20; i++) {
    for(int j = 0; j < 20; j++) {
      if(i == j) {
	Q[i][j] = new VirtualSubstitutionRate(aa[i] + aa[j], u);
      } else if(i < j) {
	Q[i][j] = new CategoryFloat(aa[i] + aa[j], rc);
      } else {
	Q[i][j] = Q[j][i];
      }
    }

    for(int j = 0; j < 20; j++) {
      if(i != j) {
	Q[i][i]->add_dependancy(Q[i][j]);
      }
    }

    RateVector* v = new RateVector("rv-" + aa[i], i, Q[i]);
    add_rate_vector(v);
  }
  finalize();
}

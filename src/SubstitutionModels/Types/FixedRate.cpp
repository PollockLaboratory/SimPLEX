#include "FixedRate.h"

#include <iostream>
#include <cstdlib>

#include "Environment.h"

extern Environment env;

FixedRate::FixedRate() {
}

void FixedRate::Initialize(int number_of_sites, std::vector<std::string> states) {
  /*
   * std::cout << "Initializing Single Probability Model" << std::endl;
   */

  float u = env.u;

  std::cout << "Protein Fixed Rate." << std::endl;

  std::array<std::string, 20> aa = {"A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"};

  std::vector<std::vector<AbstractValue*>> Q(20, std::vector<AbstractValue*>(20, NULL));

  AbstractValue* r = new FixedFloat("Rate", 0.005);
  for(int i = 0; i < 20; i++) {
    for(int j = 0; j < 20; j++) {
      if(i == j) {
	Q[i][j] = new VirtualSubstitutionRate(aa[i] + aa[j], u);
      } else {
	Q[i][j] = r;
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

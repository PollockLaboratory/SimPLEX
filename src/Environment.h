#ifndef Environment_h_
#define Environment_h_

#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <iostream>

#include "cpptoml/cpptoml.h"

using std::string;

class Environment {
public:
  Environment();
  void PrintOptions();
  void ReadOptions(int argc, char* argv[]);         // read in default and control files

  // Key values.
  int seed;                   // random number generator seed
  bool debug;                 // turns debugging on or off
  bool ancestral_sequences;   // true if the ancestral_sequences have already been sampled.
  int n; 			    // sequence length.
  int num_states;		    // number of states, not including gap.
  std::map<std::string, int> state_to_integer;

  template<typename T>
  T get(std::string);
private:
  void InitializeRandomNumberGeneratorSeed();
  void ReadTOMLfile(std::string);
  std::shared_ptr<cpptoml::table> config;
};

template<typename T>
T Environment::get(std::string option) {
  auto val = config->get_qualified_as<T>(option);
  if(val) {
    return(*val);
  } else {
    std::cerr << "Error: Option \"" << option << "\" was either not set or not the correct type." << std::endl;
    exit(EXIT_FAILURE);
  }
}

#endif

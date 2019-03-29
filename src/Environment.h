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

  template<typename T>
  std::vector<T> get_array(std::string);
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

template<typename T>
std::vector<T> Environment::get_array(std::string option) {
 std::vector<T> ret = {};
 auto nested = config;
  while(option.find(".") != std::string::npos) {
    nested = nested->get_table(option.substr(0, option.find(".")));
    option = option.substr(option.find(".") + 1);
  }
 auto vals = nested->get_array_of<T>(option);
 for(auto val : *vals) {
   ret.push_back(val);
 }
 return(ret);
}

#endif

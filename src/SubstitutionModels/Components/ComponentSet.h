#ifndef ParameterSet_h_
#define ParameterSet_h_

#include <list>
#include <vector>
#include <map>
#include <string>
// #include <iostream>
#include <fstream>

#include "AbstractComponent.h"
#include "RateVector.h"
#include "../SubstitutionModelParser.h"

struct States {
  int n; // Careful about indels.
  std::set<std::string> possible;
  std::map<std::string, int> state_to_int;
  std::map<int, std::string> int_to_state;
};

class ComponentSet {
 public:
  ComponentSet();
  void Initialize();
  void add_parameter(AbstractComponent* param);
  void add_rate_vector(RateVector* v);

  AbstractValue* realize_component(IO::raw_param);
  void create_parameters(std::list<IO::raw_param>);
  RateVector* create_rate_vector(struct States, IO::raw_rate_vector, UniformizationConstant*);
  
  bool sample();
  void accept();
  void reject();

  std::list<AbstractComponent*> get_current_parameters();
  AbstractComponent* get_current_parameter();

  void print();
  void print_dependencies();
  void refresh_dependencies();

  double get(const std::string &name);
  int size();

  void saveToFile(int gen, double l);
  std::map<AbstractComponent*, std::list<AbstractComponent*>> value_to_dependents; // Maps AbstractComponent to AbstractComponent that depend on them.
 private:
  void stepToNextParameter();

  std::list<SampleableValue*> samplable_parameters_list;
  std::list<AbstractComponent*> all_parameters_list;
  std::list<SampleableValue*>::iterator current_parameter; //Tracks the current parameter to be sampled, via an iterator across the parameter_list.

  // Dependancies.
  std::list<AbstractComponent*> get_dependent_parameters(AbstractComponent* v);

  void refreshDependancies(AbstractComponent*);

  std::map<std::string, SampleableValue*> name_to_address; //A map from the name of a parameter to the pointer of the parameter class.
  std::map<int, AbstractValue*> id_to_address; // This holds tmp original rates taken from raw_sm.
  static std::ofstream out_file;
};

#endif
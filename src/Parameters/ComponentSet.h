#ifndef ParameterSet_h_
#define ParameterSet_h_

#include <list>
#include <vector>
#include <map>
#include <string>
// #include <iostream>
#include <fstream>

#include "AbstractValue.h"
#include "RateVector.h"

class ComponentSet {
 public:
  ComponentSet();
  void Initialize();
  void add_parameter(SampleableValue* param);
  void add_rate_vector(RateVector* v);

  bool sample();
  void accept();
  void reject();

  std::list<AbstractComponent*> get_current_parameters();

  void print();
  double get(const std::string &name);
  int size();

  void saveToFile(int gen, double l);
 private:
  void stepToNextParameter();

  std::list<SampleableValue*> samplable_parameters_list;
  std::list<AbstractComponent*> all_parameters_list;
  std::list<SampleableValue*>::iterator current_parameter; //Tracks the current parameter to be sampled, via an iterator across the parameter_list.

  // Dependancies.
  std::map<AbstractComponent*, std::list<AbstractComponent*>> value_to_dependents; // Maps AbstractValues to AbstractDependentParameters that depend on them.
  std::list<AbstractComponent*> get_dependent_parameters(AbstractComponent* v);

  void setupDependancies();
  void refreshDependancies(AbstractComponent*);

  std::map<std::string, SampleableValue*> name_to_address; //A map from the name of a parameter to the pointer of the parameter class.
  static std::ofstream out_file;
};

class RateCategories : AbstractComponent {

};

#endif

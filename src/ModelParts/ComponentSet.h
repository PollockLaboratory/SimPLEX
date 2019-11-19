#ifndef ComponentSet_h_
#define ComponentSet_h_

#include <list>
#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <set>

#include "AbstractComponent.h"

/*
 * COMPONENT SET
 * This class holds a list of all components that can change through the course of the MCMC.
 * This includes:
 *   - Rate parameters.
 *   - Rate parameters that are not themselves directly sampled by are dependent on sampled rate parameters. For example virtual substitution rates.
 *   - The tree
 *   - Uniformization constant.
 */

struct SampleCounter {
  SampleableComponent* ptr; // Pointer to component.
  unsigned int freq; // Maximum sample frequency.
  unsigned int last_sample; // The last time this component was sampled
};

class ComponentSet {
public:
  ComponentSet();
  void Initialize();
  void add_parameter(AbstractComponent* param);
  void add_parameter(AbstractComponent* param, int i);

  sample_status sample();
  void accept();
  void reject();

  SampleableComponent* get_current_parameter();

  void print();
  void print_dependencies();
  void refresh_dependencies();

  double get(const std::string &name);
  int size();

  void saveToFile(int gen, double l);
  std::map<AbstractComponent*, std::list<AbstractComponent*>> value_to_dependents; // Maps AbstractComponent to AbstractComponent that depend on them.
private:
  unsigned int steps; // Number of times sample() has been called. Used as reference for differant frequencies of component samples.
  void stepToNextParameter();

  std::list<AbstractComponent*> all_parameters_list;
  std::list<SampleCounter> sampleable_parameter_list;

  std::list<SampleCounter>::iterator current_parameter; //Tracks the current parameter to be sampled, via an iterator across the parameter_list.

  // Dependancies.
  std::list<AbstractComponent*> get_dependent_parameters(AbstractComponent* v);

  void refreshDependancies(AbstractComponent*);

  std::map<std::string, Valuable*> name_to_address; //A map from the name of a parameter to the pointer of the parameter object.
  static std::ofstream out_file;
};

#endif

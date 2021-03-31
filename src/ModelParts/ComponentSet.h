#ifndef ComponentSet_h_
#define ComponentSet_h_

#include <list>
#include <map>
#include <string>

#include "AbstractComponent.h"
#include "../SubstitutionCounts.h"

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
private:
  std::list<SampleCounter> sampleable_parameter_list;
  std::map<int, AbstractComponent*> all_parameters;
  std::map<int, SubstitutionCounts*> subs_by_parameter;
  std::list<int> state_parameters;
 
  // Tracking the current parameter.
  std::list<SampleCounter>::iterator current_parameter; //Tracks the current parameter to be sampled, via an iterator across the parameter_list.
  unsigned int steps; // Number of times sample() has been called. Used as reference for differant frequencies of component samples.
  void stepToNextParameter();

  // Dependancies.
  void refresh_dependancies(AbstractComponent*);

  // Ptr to counts.
  SubstitutionCounts* counts;
public:
  ComponentSet();
  void Initialize();
  void add_parameter(AbstractComponent* param);
  void add_parameter(AbstractComponent* param, unsigned int max_sample_freq);
  void add_state_parameter(AbstractComponent* param, unsigned int max_sample_freq);
  void set_counts(SubstitutionCounts* counts);

  // Sampling.
  sample_status sample();
  void accept();
  void reject();

  // Getters.
  SampleableComponent* get_current_parameter();
  void reset_dependencies();

  // Output
  void print();
  void print_dependencies();
  void save_to_file(int gen, double l);
};

#endif

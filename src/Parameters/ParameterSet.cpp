#include "ParameterSet.h"
#include "RateVector.h"
#include "Environment.h"
#include "IO.h"

extern Environment env;
extern IO::Files files;

std::ofstream ParameterSet::out_file;

// Constructors.
ParameterSet::ParameterSet() {
  /*
   * The default parameter set constructor.
   */
}

// Setup.
void ParameterSet::Initialize() {
  current_parameter = samplable_parameters_list.begin();

  // Set up deps.
  setupDependancies();

  // Refreshes all dependancies.
  for(auto p = all_parameters_list.begin(); p != all_parameters_list.end(); ++p) {
    (*p)->refresh();
    refreshDependancies(*p);
  }

  files.add_file("parameters", env.get("parameters_out_file"), IOtype::OUTPUT);
  out_file = files.get_ofstream("parameters");

  out_file << "I,GEN,LogL";
  for(auto it = samplable_parameters_list.begin(); it != samplable_parameters_list.end(); ++it) {
    out_file << "," << (*it)->name;
  }
  for(auto it = all_parameters_list.begin(); it != all_parameters_list.end(); ++it) {
    out_file << "," << (*it)->name;
  }
  out_file << std::endl;
}

void ParameterSet::add_parameter(AbstractParameter* param) {
  /* 
   * Adds the pointer to an actual parameter onto the parameter_list, and to the
   * name_to_adress map.
   */
  samplable_parameters_list.push_back(param);

  std::string name = param->name;
  name_to_address.insert(std::make_pair(name, param));
}

void ParameterSet::add_rate_vector(RateVector* v) {
  std::vector<AbstractValue*> r = v->rates;
  for(std::vector<AbstractValue*>::iterator it = r.begin(); it != r.end(); ++it) {
    // Check if Parameter has already been seen.
    if(value_to_dependents.find(*it) == value_to_dependents.end()) {
      AbstractParameter* p = dynamic_cast<AbstractParameter*> (*it);
      if(p != NULL) {
	add_parameter(p);
      }
      value_to_dependents[*it] = {};
      all_parameters_list.push_back(*it);
    }
  }
}

void ParameterSet::setupDependancies() {
  for(auto p = all_parameters_list.begin(); p != all_parameters_list.end(); ++p) {
    std::list<AbstractValue*> deps = (*p)->get_dependancies();
    for(auto d = deps.begin(); d != deps.end(); ++d) {
      value_to_dependents[*d].push_back(*p);
    }
  }
}

void ParameterSet::refreshDependancies(AbstractValue* v) {
  std::list<AbstractValue*> deps = value_to_dependents[v];
  for(auto d = deps.begin(); d != deps.end(); ++d) {
    (*d)->refresh();
    refreshDependancies(*d);
  }
}

// Sampling.
bool ParameterSet::sample() {
  /*
   * Sample the current parameters.
   */

  // If no samplable parameters.
  if(samplable_parameters_list.empty()) {
    // Returning false will skip Metropolis Hastings step.
    return(false);
  }
  
  bool sampleType = (*current_parameter)->sample();

  try {
    refreshDependancies(*current_parameter);
  }

  catch(OutOfBoundsException &exception) {
    // A dependent parameter is out of bounds, undo the change
    // and try again.
    (*current_parameter)->undo();
    refreshDependancies(*current_parameter);
    sampleType = sample();
  }

  (*current_parameter)->refresh_host_vectors();

  return (sampleType);
}

inline void ParameterSet::stepToNextParameter() {
  /*
   * Sets the current_parameter iterator to the next sample.
   */
  ++current_parameter;
  if(current_parameter == samplable_parameters_list.end()) {
    current_parameter = samplable_parameters_list.begin();
  }
}

void ParameterSet::accept() {
  (*current_parameter)->fix();
  stepToNextParameter();
}

void ParameterSet::reject() {
  (*current_parameter)->undo();
  refreshDependancies(*current_parameter);
  (*current_parameter)->refresh_host_vectors();
  stepToNextParameter();
}

std::list<AbstractValue*> ParameterSet::get_dependent_parameters(AbstractValue* v) {
  std::list<AbstractValue*> l = {};

  std::list<AbstractValue*> deps = value_to_dependents[v];
  for(auto d = deps.begin(); d != deps.end(); ++d) {
    l.push_back(*d);
    l.splice(l.end(), get_dependent_parameters(*d));
  }
  return(l);
}

std::list<AbstractValue*> ParameterSet::get_current_parameters() {
  /*
   * Given the position of the current_parameter iterator, will return all the dependent parameters.
   * This reflects all the parameters that have changed with current sampling.
   */

  std::list<AbstractValue*> l = {*current_parameter};
  std::list<AbstractValue*> deps = get_dependent_parameters(*current_parameter);
  for(auto it = deps.begin(); it != deps.end(); ++it) {
    l.push_back(*it);
  }
  return(l);
}

void ParameterSet::print() {
  /*
   * Prints a short description of the state of the parameter_list.
   */
  std::cout << "Parameter Set - size: " << all_parameters_list.size() << std::endl;
  for(auto iter = all_parameters_list.begin(); iter != all_parameters_list.end(); ++iter) {
    (*iter)->printValue();
    std::cout << "Host vectors: ";
    for(auto i = (*iter)->host_vectors.begin(); i != (*iter)->host_vectors.end(); ++i) {
      std::cout << (*i)->name << " ";
    }
    std::cout << std::endl;
  }
}

double ParameterSet::get(const std::string &name) {
  /*
   * Will retreive the value of a parameter from the parameter set.
   */
  return name_to_address[name]->getValue();
}

int ParameterSet::size() {
  return(samplable_parameters_list.size());
}

void ParameterSet::saveToFile(int gen, double l) {
  /*
   * Saves the current parameter values to the output csv file, contained
   * in the out_file.
   */
  static int i = -1;
  ++i;
  out_file << i << "," << gen << "," << l;

  for(auto it = samplable_parameters_list.begin(); it != samplable_parameters_list.end(); ++it) {
    out_file << "," << (*it)->getValue();
  }

  for(auto it = all_parameters_list.begin(); it != all_parameters_list.end(); ++it) {
    out_file << "," << (*it)->getValue();
  }
  out_file << std::endl;
}

#include "ComponentSet.h"
#include "RateVector.h"
#include "../../Environment.h"
#include "../../IO.h"

extern Environment env;
extern IO::Files files;

std::ofstream ComponentSet::out_file;

// Constructors.
ComponentSet::ComponentSet() {
  /*
   * The default parameter set constructor.
   */
}

// Setup.

void ComponentSet::Initialize() {
  current_parameter = samplable_parameters_list.begin();

  // Refreshes all dependancies.
  for(auto p = all_parameters_list.begin(); p != all_parameters_list.end(); ++p) {
    (*p)->refresh(); // This will need some exception handeling - Virtual Sub Rate OutOfBounds etc.
    refreshDependancies(*p);
  }

  files.add_file("parameters", env.get<std::string>("OUTPUT.parameters_out_file"), IOtype::OUTPUT);
  out_file = files.get_ofstream("parameters");

  out_file << "I,GEN,LogL";
  for(auto it = samplable_parameters_list.begin(); it != samplable_parameters_list.end(); ++it) {
    out_file << "," << (*it)->get_name();
  }
  for(auto it = all_parameters_list.begin(); it != all_parameters_list.end(); ++it) {
    out_file << "," << (*it)->get_name();
  }
  out_file << std::endl;
}

void ComponentSet::add_parameter(AbstractComponent* param) {
  /* 
   * Adds the pointer to an actual parameter onto the parameter_list, and to the
   * name_to_adress map.
   */
  if(value_to_dependents.find(param) == value_to_dependents.end()) {
    // Check if Parameter has already been seen.
      SampleableValue* p = dynamic_cast<SampleableValue*> (param);
      if(p != NULL) {
	samplable_parameters_list.push_back(p);

	std::string name = p->get_name();
	name_to_address.insert(std::make_pair(name, p));
      }

      value_to_dependents[param] = {};
      all_parameters_list.push_back(param);

      // Add all deps
      std::list<AbstractComponent*> deps = param->get_dependancies();
      for(auto d = deps.begin(); d != deps.end(); ++d) {
	add_parameter(*d);
	value_to_dependents[*d].push_back(param);
      }
  }
}

void ComponentSet::add_rate_vector(RateVector* v) {
  std::vector<AbstractValue*> r = v->rates;
  for(std::vector<AbstractValue*>::iterator it = r.begin(); it != r.end(); ++it) {
    add_parameter(*it);
  }
}


AbstractComponent* ComponentSet::get_current_parameter() {
  return(*current_parameter);
}

void ComponentSet::refreshDependancies(AbstractComponent* v) {
  std::list<AbstractComponent*> deps = value_to_dependents[v];
  for(auto d = deps.begin(); d != deps.end(); ++d) {
    (*d)->refresh();
    refreshDependancies(*d);
  }
}

// Sampling.
bool ComponentSet::sample() {
  /*
   * Sample the current parameters.
   */

  bool sampleType = (*current_parameter)->sample();

  try {
    refreshDependancies(*current_parameter);
  }

  catch(OutOfBoundsException &exception) {
    // A dependent parameter is out of bounds, undo the change
    // and try again.
    (*current_parameter)->undo();
    refreshDependancies(*current_parameter);

    stepToNextParameter();
    sampleType = sample();
  }

  return (sampleType);
}

inline void ComponentSet::stepToNextParameter() {
  /*
   * Sets the current_parameter iterator to the next sample.
   */
  ++current_parameter;
  if(current_parameter == samplable_parameters_list.end()) {
    current_parameter = samplable_parameters_list.begin();
  }
}

void ComponentSet::accept() {
  (*current_parameter)->fix();
  stepToNextParameter();
}

void ComponentSet::reject() {
  (*current_parameter)->undo();
  refreshDependancies(*current_parameter);
  // (*current_parameter)->refresh_host_vectors();
  stepToNextParameter();
}

std::list<AbstractComponent*> ComponentSet::get_dependent_parameters(AbstractComponent* v) {
  std::list<AbstractComponent*> l = {};
  if(value_to_dependents[v].empty()) {
    return(std::move(l));
  } else {
    for(auto d = value_to_dependents[v].begin(); d != value_to_dependents[v].end(); ++d) {
      l.push_back(*d);
      l.splice(l.end(), get_dependent_parameters(*d));
    }
    return(l);
  }
}

std::list<AbstractComponent*> ComponentSet::get_current_parameters() {
  /*
   * Given the position of the current_parameter iterator, will return all the dependent parameters.
   * This reflects all the parameters that have changed with current sampling.
   */

  std::list<AbstractComponent*> l = {*current_parameter};
  std::list<AbstractComponent*> deps = get_dependent_parameters(*current_parameter);
  for(auto it = deps.begin(); it != deps.end(); ++it) {
    l.push_back(*it);
  }
  return(l);
}

void ComponentSet::print() {
  /*
   * Prints a short description of the state of the parameter_list.
   */
  std::cout << "Parameter Set - size: " << all_parameters_list.size() << std::endl;
  for(auto iter = all_parameters_list.begin(); iter != all_parameters_list.end(); ++iter) {
    (*iter)->print();
  }
}

void ComponentSet::print_dependencies() {
  std::cout << "Printing Dependencies." << std::endl;
  for(auto it = value_to_dependents.begin(); it != value_to_dependents.end(); ++it) {
    std::cout << (*it).first->get_name() << " [ ";
    for(auto jt = (*it).second.begin(); jt != (*it).second.end(); ++jt) {
      std::cout << (*jt)->get_name() << " ";
    }
    std::cout << "]" << std::endl;
  }

  std::cout << "SampleableParameters: " << std::endl;
  for(auto it = samplable_parameters_list.begin(); it != samplable_parameters_list.end(); ++it) {
    std::cout << (*it)->get_name() << " ";
  }
  std::cout << std::endl;
}

void ComponentSet::refresh_dependencies() {
  for(auto p = all_parameters_list.begin(); p != all_parameters_list.end(); ++p) {
    (*p)->refresh(); // This will need some exception handeling - Virtual Sub Rate OutOfBounds etc.
    refreshDependancies(*p);
  }
}

double ComponentSet::get(const std::string &name) {
  /*
   * Will retreive the value of a parameter from the parameter set.
   */
  return name_to_address[name]->getValue();
}

int ComponentSet::size() {
  return(samplable_parameters_list.size());
}

void ComponentSet::saveToFile(int gen, double l) {
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
    out_file << "," << 0.0; //(*it)->getValue();
  }
  out_file << std::endl;
}

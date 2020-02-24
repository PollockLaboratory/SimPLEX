#include "ComponentSet.h"
#include "SubstitutionModels/RateVector.h"
#include "../Environment.h"
#include "../IO/Files.h"

#include <iostream>
#include <sstream>

extern Environment env;
extern IO::Files files;

// Constructors.
ComponentSet::ComponentSet() {
  /*
   * The default parameter set constructor.
   */
  steps = 0;
}

// Setup.

void ComponentSet::Initialize() {
  current_parameter = sampleable_parameter_list.begin();

  // Reverse dependancies to set up dependents.
  for(auto p = all_parameters.begin(); p != all_parameters.end(); ++p) {
    std::list<AbstractComponent*> deps = (*p)->get_dependancies();
    for(auto d = deps.begin(); d != deps.end(); ++d) {
      (*d)->add_dependent(*p);
    }
  }

  // Setup refresh list.
  for(auto p = all_parameters.begin(); p != all_parameters.end(); ++p) {
    (*p)->setup_refresh_list();
  }

  // Refreshes all dependancies.
  refresh_all_dependencies();

  files.add_file("parameters_out", env.get<std::string>("OUTPUT.parameters_out_file"), IOtype::OUTPUT);

  std::ostringstream buffer;
  buffer << "I,GEN,LogL";
  for(auto it = all_parameters.begin(); it != all_parameters.end(); ++it) {
    buffer << "," << (*it)->get_name();
  }
  buffer << std::endl;

  files.write_to_file("parameters_out", buffer.str());
}

// Setting up.
void ComponentSet::add_parameter(AbstractComponent* param) {
  add_parameter(param, 0);
}

void ComponentSet::add_parameter(AbstractComponent* param, int max_sample_freq) {
  /* 
   * Adds the pointer to an actual parameter onto the parameter_list, and to the
   * name_to_adress map.
   */
  if(all_parameters.find(param) == all_parameters.end()) {
    // Check if Parameter has already been seen.
    SampleableComponent* p = dynamic_cast<SampleableComponent*> (param);
    if(p != NULL) {
      sampleable_parameter_list.push_back({p, max_sample_freq, 0});
    }

    all_parameters.insert(param);

    }
}

// Utils.

SampleableComponent* ComponentSet::get_current_parameter() {
  return((*current_parameter).ptr);
}

void ComponentSet::refresh_dependancies(AbstractComponent* v) {
  for(auto c = v->get_refresh_list().begin(); c != v->get_refresh_list().end(); ++c) {
    (*c)->refresh();
  }
}

void ComponentSet::refresh_all_dependencies() {
  for(auto p = all_parameters.begin(); p != all_parameters.end(); ++p) {
    refresh_dependancies(*p);
  }
}

// Sampling.
sample_status ComponentSet::sample() {
  /*
   * Sample the current parameters.
   */

  steps++;

  SampleableComponent* param = get_current_parameter();

  //std::cout << std::endl << "Start: " << param->get_name() << std::endl;
  sample_status s = param->sample();

  try {
    refresh_dependancies(param);
  }

  catch(OutOfBoundsException &exception) {
    // A dependent parameter is out of bounds, undo the change
    // and try again.
    param->undo();

    refresh_dependancies(param);

    stepToNextParameter();
    s = sample();
  }

  (*current_parameter).last_sample = steps;

  return (s);
}

inline void ComponentSet::stepToNextParameter() {
  /*
   * Sets the current_parameter iterator to the next sample.
   */
  ++current_parameter;
  if(current_parameter == sampleable_parameter_list.end()) {
    current_parameter = sampleable_parameter_list.begin();
  }

  if((steps - (*current_parameter).last_sample) + 1 < (*current_parameter).freq) {
    stepToNextParameter();
  }
}

void ComponentSet::accept() {
  get_current_parameter()->fix();
  stepToNextParameter();
}

void ComponentSet::reject() {
  get_current_parameter()->undo();
  refresh_dependancies(get_current_parameter());
  stepToNextParameter();
}

void ComponentSet::print() {
  /*
   * Prints a short description of the state of the parameter_list.
   */
  std::cout << "Component Set - size: " << all_parameters.size() << std::endl;
  for(auto iter = all_parameters.begin(); iter != all_parameters.end(); ++iter) {
    std::string sampleable;
    if(dynamic_cast<SampleableComponent*>(*iter) != nullptr) {
      sampleable = "Sampled";
    } else {
      sampleable = "Non-sampled";
    }
    std::cout << "[" << (*iter)->get_ID() << "] " << sampleable << "\t";
    (*iter)->print();
  }
}

void ComponentSet::print_dependencies() {
  //std::cout << "Printing Dependencies." << std::endl;
  //for(auto it = value_to_dependents.begin(); it != value_to_dependents.end(); ++it) {
  //  std::cout << (*it).first->get_name() << " [ ";
  //  for(auto jt = (*it).second.begin(); jt != (*it).second.end(); ++jt) {
  //    std::cout << (*jt)->get_name() << " ";
  //  }
  //  std::cout << "]" << std::endl;
  //}

  std::cout << "SampleableParameters: " << std::endl;
  for(auto it = sampleable_parameter_list.begin(); it != sampleable_parameter_list.end(); ++it) {
    std::cout << (*it).ptr->get_name() << " ";
  }
  std::cout << std::endl;
}

void ComponentSet::saveToFile(int gen, double l) {
  /*
   * Saves the current parameter values to the output csv file, contained
   * in the out_file.
   */
  static int i = -1;
  ++i;

  std::string line = std::to_string(i) + "," + std::to_string(gen) + "," + std::to_string(l);

  for(auto param = all_parameters.begin(); param != all_parameters.end(); ++param) {
    line += "," + std::to_string((*param)->record_state(gen, l));
  }

  files.write_to_file("parameters_out", line + "\n");
}

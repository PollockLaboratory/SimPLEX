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
	current_parameter = parameter_list.begin();

	// Set up deps.
	setupDependancies();
	for(auto p = parameter_list.begin(); p != parameter_list.end(); ++p) {
		refreshDependancies(*p);
	}

	files.add_file("parameters", env.get("parameters_out_file"), IOtype::OUTPUT);
	out_file = files.get_ofstream("parameters");

	out_file << "I,GEN,LogL";
	for(auto it = parameter_list.begin(); it != parameter_list.end(); ++it) {
		out_file << "," << (*it)->name;
	}
	for(auto it = dependent_parameter_list.begin(); it != dependent_parameter_list.end(); ++it) {
		out_file << "," << (*it)->name;
	}
	out_file << std::endl;
}

void ParameterSet::add_parameter(AbstractParameter* param) {
	/* 
	 * Adds the pointer to an actual parameter onto the parameter_list, and to the
	 * name_to_adress map.
	 */
	parameter_list.push_back(param);

	std::string name = param->name;
	name_to_address.insert(std::make_pair(name, param));
}

void ParameterSet::add_rate_vector(RateVector* v) {
	std::vector<AbstractValue*> r = v->rates;	
	for(std::vector<AbstractValue*>::iterator it = r.begin(); it != r.end(); ++it) {
		if(value_to_dependents.find(*it) == value_to_dependents.end()) {
			AbstractParameter* p = dynamic_cast<AbstractParameter*> (*it);
			if(p != NULL) {
				value_to_dependents[p] = {};
				add_parameter(p);
			} else {
				AbstractDependentParameter* hp = dynamic_cast<AbstractDependentParameter*> (*it);
				value_to_dependents[hp] = {};
				dependent_parameter_list.push_back(hp);
			}
		}
	}
}

void ParameterSet::setupDependancies() {
	for(auto p = dependent_parameter_list.begin(); p != dependent_parameter_list.end(); ++p) {
		std::vector<AbstractValue*> deps = (*p)->get_dependancies();	
		for(auto d = deps.begin(); d != deps.end(); ++d) {
			value_to_dependents[(*d)].push_back(*p);
		}
	}
}

void ParameterSet::refreshDependancies(AbstractValue* v) {
	std::list<AbstractDependentParameter*> deps = value_to_dependents[v];
	for(auto d = deps.begin(); d != deps.end(); ++d) {
		(*d)->refresh();
		refreshDependancies(*d);
	}
}

// Sampling.
bool ParameterSet::sample() {
	/*
	 * Will sample the current parameters.
	 */
	bool sampleType = (*current_parameter)->sample();
	refreshDependancies(*current_parameter);
	return (sampleType);
}

inline void ParameterSet::stepToNextParameter() {
	/*
	 * Sets the current_parameter iterator to the next sample.
	 */
	++current_parameter;
	if(current_parameter == parameter_list.end()) {
		current_parameter = parameter_list.begin();
	}
}

void ParameterSet::accept() {
	(*current_parameter)->fix();
	stepToNextParameter();
}

void ParameterSet::reject() {
	(*current_parameter)->undo();
	refreshDependancies(*current_parameter);
	stepToNextParameter();
}

void ParameterSet::print() {
	/*
	 * Prints a short description of the state of the parameter_list.
	 */
	std::cout << "Parameter Set - size: " << parameter_list.size() << std::endl;
	for(auto iter = parameter_list.begin(); iter != parameter_list.end(); ++iter) {
		(*iter)->printValue();
	}
	std::cout << "Dependent Parameter Set - size: " << dependent_parameter_list.size() << std::endl;
	for(auto iter = dependent_parameter_list.begin(); iter != dependent_parameter_list.end(); ++iter) {
		(*iter)->printValue();
	}
}

double ParameterSet::get(const std::string &name) {
	/*
	 * Will retreive the value of a parameter from the parameter set.
	 */
	return name_to_address[name]->getValue();
}

int ParameterSet::size() {
	return(parameter_list.size());
}

void ParameterSet::saveToFile(int gen, double l) {
	/*
	 * Saves the current parameter values to the output csv file, contained
	 * in the out_file.
	 */
	static int i = -1;
	++i;
	out_file << i << "," << gen << "," << l;
	for(auto it = parameter_list.begin(); it != parameter_list.end(); ++it) {
		out_file << "," << (*it)->getValue();
	}
	for(auto it = dependent_parameter_list.begin(); it != dependent_parameter_list.end(); ++it) {
		out_file << "," << (*it)->getValue();
	}
	out_file << std::endl;
}

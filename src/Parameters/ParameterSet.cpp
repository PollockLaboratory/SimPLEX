#include "ParameterSet.h"
#include "RateVector.h"

// Constructors.
ParameterSet::ParameterSet() {
	/*
	 * The default parameter set constructor.
	 */
}

// Setup.
void ParameterSet::Initialize(std::ofstream* &out_file_buffer) {
	current_parameter = parameter_list.begin();
	out_stream_buffer = out_file_buffer;

	// Set up deps.
	setupDependancies();
	for(auto p = parameter_list.begin(); p != parameter_list.end(); ++p) {
		refreshDependancies(*p);
	}

	AddHeaderToFile();
	RecordStateToFile();
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
void ParameterSet::sample() {
	/*
	 * Will sample the current parameters.
	 */
	(*current_parameter)->sample();
	refreshDependancies(*current_parameter);
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

void ParameterSet::AddHeaderToFile() {
	/*
	 * Adds the columns names to the output csv file.
	 */
	std::list<AbstractParameter*>::iterator iter = parameter_list.begin();
	while(iter != parameter_list.end()) {
		if(iter != parameter_list.begin()) {
			*out_stream_buffer << ", ";

		}
		*out_stream_buffer << (*iter)->name;
		++iter;
	}
	*out_stream_buffer << std::endl;
}

void ParameterSet::RecordStateToFile(){
	/*
	 * Saves the current parameter values to the output csv file, contained
	 * in the out_stream_buffer.
	 */
	std::list<AbstractParameter*>::iterator iter = parameter_list.begin();
	while(iter != parameter_list.end()) {
		if(iter != parameter_list.begin()) {
			*out_stream_buffer << ", ";

		}
		*out_stream_buffer << (*iter)->getValue();
		++iter;
	}
	*out_stream_buffer << std::endl;
}

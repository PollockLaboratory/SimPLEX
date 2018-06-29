#include "ParameterSet.h"

ParameterSet::ParameterSet() {
	/*
	 * The default parameter set constructor.
	 */
}

void ParameterSet::Initialize(std::ofstream* &out_file_buffer) {
	current_parameter = parameter_list.begin();
	out_stream_buffer = out_file_buffer;

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

void ParameterSet::sample() {
	/*
	 * Will sample the current parameters.
	 */
	(*current_parameter)->sample();
	
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
	stepToNextParameter();
}

void ParameterSet::print() {
	/*
	 * Prints a short description of the state of the parameter_list.
	 */
	std::cout << "Parameter Set - size: " << parameter_list.size() << std::endl;
	std::list<AbstractParameter*>::iterator iter = parameter_list.begin();
	while(iter != parameter_list.end()) {
		(*iter)->printValue();
		++iter;
	}
}

double ParameterSet::get(const std::string &name) {
	/*
	 * Will retreive the value of a parameter from the parameter set.
	 */
	return name_to_address[name]->getValue();
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





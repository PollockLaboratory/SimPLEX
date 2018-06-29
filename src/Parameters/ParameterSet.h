#ifndef ParameterSet_h_
#define ParameterSet_h_

#include <list>
#include <map>
#include <string>
#include <iostream>
#include <fstream>

#include "AbstractParameter.h"

class ParameterSet {
	public:
		ParameterSet();
		void Initialize(std::ofstream* &out_file_buffer);
		void add_parameter(AbstractParameter*);

		void sample();
		void accept();
		void reject();

		void print();
		double get(const std::string &name);
		void RecordStateToFile();
	private:
		void stepToNextParameter();
		void AddHeaderToFile();

		std::list<AbstractParameter*> parameter_list;
		std::list<AbstractParameter*>::iterator current_parameter; //Tracks the current parameter to be sampled, via an iterator across the parameter_list.
		std::map<std::string, AbstractParameter*> name_to_address; //A map from the name of a parameter to the pointer of the parameter class.
		std::ofstream* out_stream_buffer;
};

#endif

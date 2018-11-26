#ifndef ParameterSet_h_
#define ParameterSet_h_

#include <list>
#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <fstream>

#include "AbstractValue.h"
#include "RateVector.h"

class ParameterSet {
	public:
		ParameterSet();
		void Initialize();
		void add_parameter(AbstractParameter* param);
		void add_rate_vector(RateVector* v);

		bool sample();
		void accept();
		void reject();

		std::list<AbstractValue*> get_current_parameters();

		void print();
		double get(const std::string &name);
		int size();

		void saveToFile(int gen, double l);
	private:
		void stepToNextParameter();

		std::list<AbstractParameter*> parameter_list;
		std::list<AbstractDependentParameter*> dependent_parameter_list;
		std::list<AbstractParameter*>::iterator current_parameter; //Tracks the current parameter to be sampled, via an iterator across the parameter_list.
	
		// Dependancies.
		std::map<AbstractValue*, std::list<AbstractDependentParameter*>> value_to_dependents; // Maps AbstractValues to AbstractDependentParameters that depend on them.
		std::list<AbstractDependentParameter*> get_dependent_parameters(AbstractValue* v);

		void setupDependancies();
		void refreshDependancies(AbstractValue*);

		std::map<std::string, AbstractParameter*> name_to_address; //A map from the name of a parameter to the pointer of the parameter class.
		static std::ofstream out_file;
};

#endif

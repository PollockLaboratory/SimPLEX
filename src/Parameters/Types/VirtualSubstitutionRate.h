#ifndef VirtualSubstitutionRate_h_
#define VirtualSubstitutionRate_h_

#include "AbstractValue.h"
#include <string>

class VirtualSubstitutionRate : public AbstractDependentParameter {
	public:
		VirtualSubstitutionRate(std::string parameter_name, double unif);
		virtual double getValue(); 	
		virtual void printValue();
		virtual void refresh();
		virtual void add_dependancy(AbstractValue*);
	private:
		double value;
		double u; //Uniformization constant.
};

#endif

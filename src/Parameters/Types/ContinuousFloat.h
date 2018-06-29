#ifndef ContinuousFloat_h_
#define ContinuousFloat_h_

#include <string>
#include <iostream>

#include "../AbstractParameter.h"

class ContinuousFloat : public AbstractParameter {
	public:
		ContinuousFloat(std::string, double, double);
		virtual void printValue();
		virtual void sample();
		virtual double getValue();
		virtual void undo();
		virtual void fix();
	private:
		double value;
		double std_dev;
		double previous_value;
		
};

#endif

#ifndef AbstractParameter_h_
#define AbstractParameter_h_

#include <string>

class AbstractParameter {
	public:
		AbstractParameter(std::string);
		std::string name;

		virtual void sample() = 0;
		virtual void undo() = 0;
		virtual void fix() = 0;

		virtual double getValue() = 0;

		virtual void printValue() = 0;
	protected:
		bool fixedQ; //Indicates whether the state of this parameter is fixed or still being trialed.
};



#endif

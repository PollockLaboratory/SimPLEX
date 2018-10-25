#ifndef AbstractValue_h_
#define AbstractValue_h_

#include <string>
#include <vector>

class AbstractValue {
	public:
		AbstractValue(std::string parameter_name) {
			name = parameter_name;
		}

		std::string name;
		virtual double getValue() = 0;
		virtual void printValue() = 0;
};

class AbstractParameter : public AbstractValue {
	public:
		AbstractParameter(std::string parameter_name) : AbstractValue(parameter_name) {
			fixedQ = true;
		}

		virtual void sample() = 0;
		virtual void undo() = 0;
		virtual void fix() = 0;
	protected:
		bool fixedQ; //Indicates whether the state of this parameter is fixed or still being trialed.
};

class AbstractDependentParameter : public AbstractValue {
	public:
		AbstractDependentParameter(std::string parameter_name) : AbstractValue(parameter_name) {}
		virtual void refresh() = 0;
		virtual void add_dependancy(AbstractValue*) = 0;
		std::vector<AbstractValue*> get_dependancies() { return(dependent_values); }
	protected:
		std::vector<AbstractValue*> dependent_values;
};

#endif

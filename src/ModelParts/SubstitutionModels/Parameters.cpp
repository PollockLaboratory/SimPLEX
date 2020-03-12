#include "Parameters.h"
#include <stdlib.h> //This gives rand.
#include <limits>

// CONTINUOUS FLOAT

ContinuousFloat::ContinuousFloat(std::string name, double initial_value = 0.0, double initial_std_dev = 1.0) : SampleableValue(name), value(initial_value), std_dev(initial_std_dev) {
  /*
   * The default constructor for the Continuous Float parameter class.
   */

  previous_value = initial_value;
}

void ContinuousFloat::print() {
  std::cout << "Continuous float - " << name << ": " << value
	    << " [" << lower_bound->get_description()
	    << "->" << upper_bound->get_description() << "]" << std::endl;
}

sample_status ContinuousFloat::sample() {
  previous_value = value;
  fixedQ = false;

  double r = ((rand() % 10000) / 10000.0) - 0.5;
  value = value + (r * std_dev);

  double lb = lower_bound->get_value();
  double ub = upper_bound->get_value();

  while(value < lb or value > ub) {
    if(value < lb) {
      value = (2*lb) - value;
    }

    if(value > ub) {
      value = (2*ub) - value;
    }
  }

  return(sample_status({true, true, false}));
}

const double& ContinuousFloat::getValue() {
  return(value);
}

const double& ContinuousFloat::getOldValue() { 
  return(previous_value);
}

void ContinuousFloat::undo() {
  value = previous_value;
  fixedQ = true;
}

void ContinuousFloat::fix() {
  previous_value = value;
  fixedQ = true;
}

void ContinuousFloat::refresh() {
}

double ContinuousFloat::record_state(int gen, double l) {
  return(getValue());
}

std::string ContinuousFloat::get_type() {
  return("CONTINUOUS_FLOAT");
}

std::ostream& operator<<(std::ostream& os, const ContinuousFloat& cf) {
  os << "[ContinuousFloat-" << cf.name << "]";
  return(os);
}

// DISCRETE FLOAT
RateCategories::RateCategories(std::string name, std::vector<Valuable*> categories) : AbstractComponent(name) {
  values = categories;
  n = values.size();
  for(auto it = values.begin(); it != values.end(); ++it) {
    AbstractComponent* component = dynamic_cast<AbstractComponent*>(*it);
    if(component != nullptr) {
      add_dependancy(component);
    } else {
      std::cerr << "Error: cannot add value dependancies to RateCategories." << std::endl;  
    }
  }
}

void RateCategories::fix() {
}

void RateCategories::refresh() {
}

double RateCategories::record_state(int gen, double l) {
  return(0.0);
}

void RateCategories::print() {
  std::cout << "Rate Categories - " << name << ": [ ";
  for(int i = 0; i < n; i++) {
    AbstractComponent* component = dynamic_cast<AbstractComponent*>(values[i]);
    std::cout << component->name << " ";
  }
  std::cout << "]" << std::endl;
}

double RateCategories::operator[](int i) {
  return(values[i]->getValue());
}

std::string RateCategories::get_type() {
  return("RATE_CATEGORIES");
}

DiscreteFloat::DiscreteFloat(std::string name, RateCategories* categories) : SampleableValue(name) {
  /*
   * The default constructor for the Continuous Float parameter class.
   */
  rc = categories;
  add_dependancy(rc);
  n = rc->n;
  i = 0;
  prev_i = 0;
  value = (*rc)[i]; 
  previous_value = value;

  if(value == 0) {
    std::cerr << "Error: category value is 0." << std::endl;
    std::cerr << value << " : " << i << std::endl;
    exit(EXIT_FAILURE);
  }
}

// Utils
void DiscreteFloat::print() {
  std::cout << "Category float - " << name << ": " << value << std::endl;
}

sample_status DiscreteFloat::sample() {
  fixedQ = false;

  prev_i = i;

  if(i == 0) {
    i = 1;
  } else if(i == n-1) {
    i = i - 1;
  } else {
    int r  = rand() % 2;
    if(r == 0) {
      i = i - 1;
    } else {
      i = i + 1;
    }
  }

  value = (*rc)[i];

  if(value == 0) {
    std::cerr << "Error: category value is 0." << std::endl;
    std::cerr << value << " : " << i << std::endl;
    exit(EXIT_FAILURE);
  }
  
  return(sample_status({true, true, false}));
}

const double& DiscreteFloat::getValue() {
  return(value);
}

const double& DiscreteFloat::getOldValue() {
  return(previous_value);
}

void DiscreteFloat::undo() {
  i = prev_i;
  value = (*rc)[i];
  fixedQ = true;
}

void DiscreteFloat::fix() {
  //value = (*rc)[i];
  prev_i = i;
  previous_value = value;
  fixedQ = true;
}

void DiscreteFloat::refresh() {
    value = (*rc)[i];
}

double DiscreteFloat::record_state(int gen, double l) {
  return(getValue());
}

std::string DiscreteFloat::get_type() {
  return("DISCRETE_FLOAT");
}

// FIXED FLOAT

FixedFloat::FixedFloat(std::string parameter_name, double v) : StaticValue(parameter_name) {
  value = v;
}

const double& FixedFloat::getValue() {
  return(value);
}

const double&  FixedFloat::getOldValue() {
  return(value);
}

void FixedFloat::print() {
  std::cout << "FixedFloat - " << name << ": " << value << std::endl;
}

void FixedFloat::fix() {
}

void FixedFloat::refresh() {	
}

std::string FixedFloat::get_type() {
  return("FIXED_FLOAT");
}

// VIRTUAL SUBSTITUTION RATE

VirtualSubstitutionRate::VirtualSubstitutionRate(AbstractComponent* parameter, Valuable* unif) : StaticValue(parameter) {
  u = unif;

  UniformizationConstant* uniform = dynamic_cast<UniformizationConstant*>(u);
  if(uniform != nullptr) {
    uniform->add_VirtualSubstitutionRate(this);
    this->add_dependancy(uniform);
  }
  // This should be more throughtfully set.
  value = 0.232323;
  previous_value = value;
}

const double& VirtualSubstitutionRate::getValue() {
  return(value);
}

const double& VirtualSubstitutionRate::getOldValue() {
  return(previous_value);
}

void VirtualSubstitutionRate::print() {
  std::cout << "VirtualSubstitutionRate - " << name << ": " << value << std::endl;
}

std::string VirtualSubstitutionRate::get_type() {
  return("VIRTUAL_SUBSTITUTION_RATE");
}

void VirtualSubstitutionRate::fix() {
  previous_value = value;
}

void VirtualSubstitutionRate::refresh() {
  double total = 0.0;

  int i = 0;
  //std::cout << "rates [ " << std::endl;
  for(auto it = dependent_rates.begin(); it != dependent_rates.end(); ++it) {
    // std::cout << (*it)->getValue() << " ";
    total += (*it)->getValue();
  }
  //std::cout << "]" << std::endl;
  
  value = u->getValue() - total;

  //std::cout << previous_value << "->" << value << std::endl;
  if(value <= 0.0 || value > 1.0) {
    throw OutOfBoundsException("VirtualSubstitutionRate out of bounds.");
  }
}

void VirtualSubstitutionRate::add_rate(Valuable* v) {
  AbstractComponent* c = dynamic_cast<AbstractComponent*>(v);
  if(c == nullptr) {
    std::cerr << "Error: in VirtualSubstitutionRate::add_rate, Value is not AbstractComponent." << std::endl;
    exit(EXIT_FAILURE);
  }
  this->add_dependancy(c);
  dependent_rates.push_back(v);
}



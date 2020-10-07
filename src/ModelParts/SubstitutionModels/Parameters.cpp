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

const double& ContinuousFloat::get_value() {
  return(value);
}

const double& ContinuousFloat::get_old_value() { 
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

std::string ContinuousFloat::get_state_header() {
  return(name);
}

std::string ContinuousFloat::get_state() {
  return(std::to_string(get_value()));
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
      exit(EXIT_FAILURE);
    }
  }

  // Check that that no value is equal to 0, and they in the right order.
  float prev_val = 0.0;
  float val;
  for(auto it = values.begin(); it != values.end(); ++it) {
    val = (*it)->get_value();
    if(val <= prev_val) {
      if(val == 0.0) {
	std::cerr << "Error in RateCategories construction: initial categories cannot be set to a value of 0.0" << std::endl;
	exit(EXIT_FAILURE);
      }
      std::cerr << "Error in RateCategories construction: initial rate categories are not ordered (from lowest to highest)." << std::endl;
      exit(EXIT_FAILURE); 
    }
    prev_val = val;
  }
}

void RateCategories::fix() {
}

void RateCategories::refresh() {
}

std::string RateCategories::get_state_header() {
  return(name);
}

std::string RateCategories::get_state() {
  return("n/a");
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
  return(values[i]->get_value());
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
    std::cerr << "Error: attempted sampling to a category at a value equal to 0.0" << std::endl;
    exit(EXIT_FAILURE);
  }
  
  return(sample_status({true, true, false}));
}

const double& DiscreteFloat::get_value() {
  return(value);
}

const double& DiscreteFloat::get_old_value() {
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

std::string DiscreteFloat::get_state_header() {
  return(name + "," + name + "-category");
}

std::string DiscreteFloat::get_state() {
  return(std::to_string(get_value()) + "," + std::to_string(i));
}

std::string DiscreteFloat::get_type() {
  return("DISCRETE_FLOAT");
}

// FIXED FLOAT

FixedFloat::FixedFloat(std::string parameter_name, double v) : NonSampleableValue(parameter_name) {
  value = v;
}

const double& FixedFloat::get_value() {
  return(value);
}

const double&  FixedFloat::get_old_value() {
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

// Arithmatic

double arith_add(double v1, double v2) {
  return(v1 + v2);
};

double arith_subtract(double v1, double v2) {
  return(v1 - v2);
};

double arith_multiply(double v1, double v2) {
  return(v1 * v2);
}

double arith_divide(double v1, double v2) {
  return(v1 / v2);
}

int Arithmatic::idc = 0;

void Arithmatic::setup(Valuable* v1, Valuable* v2) {
  switch(op) {
  case ADDITION:
    arith_func = &arith_add;
    break;
  case SUBTRACTION:
    arith_func = &arith_subtract;
    break;
  case MULTIPLICATION:
    arith_func = &arith_multiply;
    break;
  case DIVISION:
    arith_func = &arith_divide;
    break;
  }
  
  AbstractComponent* c1 = dynamic_cast<AbstractComponent*>(v1);
  AbstractComponent* c2 = dynamic_cast<AbstractComponent*>(v2);

  this->add_dependancy(c1);
  this->add_dependancy(c2);

  value = arith_func(v1->get_value(), v2->get_value());
  previous_value = value;
}

Arithmatic::Arithmatic(std::string name, arith_op op, Valuable* v1, Valuable* v2) : NonSampleableValue(name), op(op), v1(v1), v2(v2) {
  setup(v1, v2);
}

Arithmatic::Arithmatic(arith_op op, Valuable* v1, Valuable* v2) : NonSampleableValue("Arithmatic-" + std::to_string(idc)), op(op), v1(v1), v2(v2) {
  idc++;
  setup(v1, v2);
}

const double& Arithmatic::get_value() {
  return(value);
}

const double& Arithmatic::get_old_value() {
  return(previous_value);
}

void Arithmatic::print() {
  std::string op_char = "";
  switch(op) {
  case ADDITION:
    op_char = "+";
    break;
  case SUBTRACTION:
    op_char = "-";
    break;
  case MULTIPLICATION:
    op_char = "*";
    break;
  case DIVISION:
    op_char = "/";
    break;
  }

  std::string name1 = dynamic_cast<AbstractComponent*>(v1)->name;
  std::string name2 = dynamic_cast<AbstractComponent*>(v2)->name;

  std::cout << "ARITH[" << name1 << op_char << name2 << "] " << name << ": " << value << std::endl;
}

void Arithmatic::fix() {
  previous_value = value;
}

void Arithmatic::refresh() {
  value = arith_func(v1->get_value(), v2->get_value());
}

std::string Arithmatic::get_type() {
  return("Arithmatic");
}

// VIRTUAL SUBSTITUTION RATE

VirtualSubstitutionRate::VirtualSubstitutionRate(std::string name) : NonSampleableValue(name) {
  u = nullptr;

  value = 0.232323;
  previous_value = value;
}

const double& VirtualSubstitutionRate::get_value() {
  return(value);
}

const double& VirtualSubstitutionRate::get_old_value() {
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

  //std::cout << "rates [ ";
  for(auto it = dependent_rates.begin(); it != dependent_rates.end(); ++it) {
    //std::cout << (*it)->get_value() << " ";
    total += (*it)->get_value();
  }
  //std::cout << "]" << std::endl;
  
  value = u->get_value() - total;

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

void VirtualSubstitutionRate::set_u(Valuable *unif) {
  u = unif;

  // If Uniformization constant is dynamic, rather that simply a fixed float.
  UniformizationConstant* uniform = dynamic_cast<UniformizationConstant*>(u);
  if(uniform != nullptr) {
    uniform->add_VirtualSubstitutionRate(this);
    this->add_dependancy(uniform);
  }
}


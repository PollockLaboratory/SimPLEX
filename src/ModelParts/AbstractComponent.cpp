#include <iostream>
#include <limits>

#include "AbstractComponent.h"
#include "SubstitutionModels/RateVector.h"
#include "../Environment.h"

extern Environment env;

// ABSTRACT COMPONENT

AbstractComponent::AbstractComponent(std::string name) : name(name) {
  static int idc = -1;
  idc++;
  ID = idc;
  refresh_list = {};
  valuable_dependents = {};
}

AbstractComponent::AbstractComponent(AbstractComponent* parameter) {
  ID = parameter->get_ID();
  name = parameter->get_name();

  refresh_list = {};
  valuable_dependents = {};
  delete parameter;
}

int AbstractComponent::get_ID() {
  return(ID);
}

std::string AbstractComponent::get_name() {
  return(name);
}

void AbstractComponent::add_dependancy(AbstractComponent* v) {
  dependencies.push_back(v);
}

const std::list<AbstractComponent*>& AbstractComponent::get_dependancies() {
  return(dependencies);
}

void AbstractComponent::add_dependent(AbstractComponent* v) {
  // Check if it already exists in the list.
  dependents.push_back(v);
}

const std::list<AbstractComponent*>& AbstractComponent::get_dependents() {
  return(dependents);
}

std::list<AbstractComponent*> AbstractComponent::next_dependents(std::list<AbstractComponent*> components, std::set<AbstractComponent*>& previous) {
  std::list<AbstractComponent*> ret = {};
  for(auto c = components.begin(); c != components.end(); ++c) {
    for(auto d = (*c)->get_dependents().begin(); d != (*c)->get_dependents().end(); ++d) {
      if(previous.find(*d) == previous.end()) {
	// There are likely duplicates in refresh_list.
	ret.push_back(*d);
	refresh_list.push_back(*d);
	Valuable* val = dynamic_cast<Valuable*>(*d);
	if(val != nullptr) {
	  valuable_dependents.push_back(val);
	}
	previous.insert(*d);
      }
    }
  }

  return(ret);
}

void AbstractComponent::setup_refresh_list() {
  // This function is naive to the way components are related to one another.
  std::list<AbstractComponent*> components = {this};
  refresh_list.push_back(this);

  Valuable* val = dynamic_cast<Valuable*>(this);
  if(val != nullptr) {
    valuable_dependents.push_back(val);
  }

  while(not components.empty()) {
    std::set<AbstractComponent*> previous = {};
    for(auto it = components.begin(); it != components.end(); ++it) {
      previous.insert(*it);
    }
    
    components = next_dependents(components, previous);
  }
};

const std::list<AbstractComponent*>& AbstractComponent::get_refresh_list() {
  return(refresh_list);
}

const std::list<Valuable*>& AbstractComponent::get_valuable_dependents() {
  return(valuable_dependents);
}

// ABSTRACT VALUES
Valuable::Valuable() {
  host_vectors = {};
}

void Valuable::add_host_vector(RateVector* rv, int pos) {
  host_vectors.push_back(rv_loc {rv, pos});
}

std::list<rv_loc> Valuable::get_host_vectors() {
  return(host_vectors);
}

void Valuable::print() {
  std::cout << this->getValue();
}

// SampleableComponent.

SampleableComponent::SampleableComponent(std::string name) : AbstractComponent(name) {
  fixedQ = true;
  dependencies = {};
}

// Constaints
double inf_d = std::numeric_limits<double>::infinity();
double neg_inf_d = -inf_d;

FixedConstraint::FixedConstraint(double value) : value(value) {
}

double FixedConstraint::get_value() const {
  return(value);
}

std::string FixedConstraint::get_description() const {
  if(value == inf_d) {
    return("inf");
  } else if(value == neg_inf_d) {
    return("-inf");
  } else {
    return(std::to_string(value));
  }
}

DynamicConstraint::DynamicConstraint(Valuable* value) : value(value) {
}

double DynamicConstraint::get_value() const {
  return(value->getValue());
}

std::string DynamicConstraint::get_description() const {
  AbstractComponent* component = dynamic_cast<AbstractComponent*>(value);
  return(component->name);
}

// SampleableValue.
FixedConstraint initial_lower_bound = FixedConstraint(neg_inf_d);
FixedConstraint initial_upper_bound = FixedConstraint(inf_d);

SampleableValue::SampleableValue(std::string name) : SampleableComponent(name), Valuable() {
  lower_bound = &initial_lower_bound;
  upper_bound = &initial_upper_bound;
}

void SampleableValue::set_lower_boundary(BaseConstraint* constraint) {
  lower_bound = constraint;
}

void SampleableValue::set_upper_boundary(BaseConstraint* constraint) {
  upper_bound = constraint;
}

// StaticValue.

StaticValue::StaticValue(std::string name) : AbstractComponent(name), Valuable() {
}

StaticValue::StaticValue(AbstractComponent* parameter) : AbstractComponent(parameter), Valuable() {
}

double StaticValue::record_state(int gen, double l) {
  return(getValue());
}

// UniformizationConstant.
// This is here because it is a special parameter.

UniformizationConstant::UniformizationConstant(double initial_value) : SampleableComponent("U"), value(initial_value), previous_value(initial_value) {
  threshold = env.get<double>("UNIFORMIZATION.threshold");
  max_step = env.get<double>("UNIFORMIZATION.max_step");
}

void UniformizationConstant::print() {
  std::cout << "DynamicUniformizationConstant: " << value << std::endl;
}

void UniformizationConstant::set_initial() {
  double min = 1.0;
  for(auto it = vsrs.begin(); it != vsrs.end(); ++it) {
    double v = (*it)->getValue();
    if(v < min) {
      min = v;
    }
  }
  value = 1.05 - min;
}

sample_status UniformizationConstant::sample() {
  previous_value = value;
  double min = 1.0;
  for(auto it = vsrs.begin(); it != vsrs.end(); ++it) {
    double v = (*it)->getValue();
    if(v < min) {
      min = v;
    }
  }

  double spare = min - threshold;
  if(spare > 0.0) {
    // Reduce uniformizationconstant.
    if(spare < max_step) {
      value = value - spare;
    } else {
      value = value - max_step;
    } 
  } else {
    value = value - spare;
  }
 
  return(sample_status({false, true, true}));
}

const double& UniformizationConstant::getValue() {
  return(value);
}

const double& UniformizationConstant::getOldValue() {
  return(previous_value);
}

void UniformizationConstant::undo() {
  value = previous_value;
}

void UniformizationConstant::fix() {
}

void UniformizationConstant::refresh() {
}

double UniformizationConstant::record_state(int gen, double l) {
  return(getValue());
}

std::string UniformizationConstant::get_type() {
  return("UNIFORMIZATION_CONSTANT");
}

void UniformizationConstant::add_VirtualSubstitutionRate(Valuable* v) {
  // These are the rates that if they get close to 0 or if they are all far from 0,
  // will change the uniformization rate.
  vsrs.push_back(v);
}

std::ostream& operator<<(std::ostream& os, const UniformizationConstant& u) {
  os << "[UniformizationConstant-" << u.value << "]";
  return(os);
}

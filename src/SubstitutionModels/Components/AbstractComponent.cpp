#include "../Environment.h"

#include "AbstractComponent.h"
#include "RateVector.h"
#include <iostream>

extern Environment env;

// ABSTRACT COMPONENT

AbstractComponent::AbstractComponent(std::string name) : name(name) {
  static int IDc = 0;
  IDc++;
  ID = IDc;
  host_vectors = {};
}

int AbstractComponent::get_ID() {
  return(ID);
}

std::string AbstractComponent::get_name() {
  return(name);
}

void AbstractComponent::add_dependancy(AbstractComponent* v) {
  dependent_values.push_back(v);
}

const std::list<AbstractComponent*>& AbstractComponent::get_dependancies() {
  return(dependent_values);
}

// ABSTRACT VALUES

AbstractValue::AbstractValue(std::string name) : AbstractComponent(name) {
  host_vectors = {};
  dependent_values = {};
}

void AbstractValue::add_host_vector(RateVector* rv, int pos) {
  host_vectors.push_back(rv_loc {rv, pos});
}

std::list<rv_loc> AbstractValue::get_host_vectors() {
  return(host_vectors);
}

// UniformizationConstant.
// This is here because it is a special parameter.

UniformizationConstant::UniformizationConstant() : SampleableValue("U"), value(1.0), previous_value(1.0) {
  threshold = env.get<double>("UNIFORMIZATION.threshold");
  max_step = env.get<double>("UNIFORMIZATION.max_step");
}

void UniformizationConstant::print() {
  std::cout << "UniformizationConstant: " << value << std::endl;
}

void UniformizationConstant::set_initial() {
  double min = 1.0;
  std::cout << "Setting initial unif: " << std::endl;
  for(auto it = vsrs.begin(); it != vsrs.end(); ++it) {
    double v = (*it)->getValue();
    std::cout << v << " ";
    if(v < min) {
      min = v;
    }
  }
  std::cout << std::endl;
  value = 1.05 - min;
}

bool UniformizationConstant::sample() {
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
 
  return(false);
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

void UniformizationConstant::add_VirtualSubstitutionRate(AbstractValue* v) {
  // These are the rates that if they get close to 0 or if they are all far from 0,
  // will change the uniformization rate.
  vsrs.push_back(v);
}

std::ostream& operator<<(std::ostream& os, const UniformizationConstant& u) {
  os << "[UniformizationConstant-" << u.value << "]";
  return(os);
}





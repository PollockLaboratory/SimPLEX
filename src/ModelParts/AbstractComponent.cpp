#include <iostream>

#include "AbstractComponent.h"
#include "SubstitutionModels/RateVector.h"
#include "../Environment.h"

extern Environment env;

// ABSTRACT COMPONENT

AbstractComponent::AbstractComponent(std::string name) : name(name) {
  static int idc = -1;
  idc++;
  ID = idc;
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
  dependent_values = {};
}

// SampleableValue.

SampleableValue::SampleableValue(std::string name) : SampleableComponent(name), Valuable() {
}

// StaticValue.

StaticValue::StaticValue(std::string name) : AbstractComponent(name), Valuable() {
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

void UniformizationConstant::add_VirtualSubstitutionRate(Valuable* v) {
  // These are the rates that if they get close to 0 or if they are all far from 0,
  // will change the uniformization rate.
  vsrs.push_back(v);
}

std::ostream& operator<<(std::ostream& os, const UniformizationConstant& u) {
  os << "[UniformizationConstant-" << u.value << "]";
  return(os);
}

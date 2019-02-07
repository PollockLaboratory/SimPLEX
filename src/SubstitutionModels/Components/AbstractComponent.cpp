#include "AbstractComponent.h"
#include "RateVector.h"
#include <iostream>

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
  host_vectors.push_back(std::pair<RateVector*, int>(rv, pos));
}

std::list<std::pair<RateVector*, int>> AbstractValue::get_host_vectors() {
  return(host_vectors);
}



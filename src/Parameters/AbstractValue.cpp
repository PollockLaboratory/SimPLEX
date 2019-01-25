#include <AbstractValue.h>
#include <RateVector.h>
#include <iostream>

// ABSTRACT COMPONENT

AbstractComponent::AbstractComponent(std::string name) : name(name) {
  static int IDc = 0;
  IDc++;
  ID = IDc;
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

std::list<AbstractComponent*> AbstractComponent::get_dependancies() {
  return(dependent_values);
}

// ABSTRACT VALUES

AbstractValue::AbstractValue(std::string name) : AbstractComponent(name) {
  host_vectors = {};
  dependent_values = {};
}

void AbstractValue::add_host_vector(RateVector* rv) {
  host_vectors.push_back(rv);
}

void AbstractValue::refresh_host_vectors() {
  for(auto it = host_vectors.begin(); it != host_vectors.end(); ++it) {
    (*it)->update_single_logLikelihood(ID);
  }
}



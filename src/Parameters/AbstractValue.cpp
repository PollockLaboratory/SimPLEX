#include <AbstractValue.h>
#include <RateVector.h>
#include <iostream>

AbstractValue::AbstractValue(std::string parameter_name) : name(parameter_name) {
  static int IDc = 0;
  IDc++;
  ID = IDc;
  host_vectors = {};
  dependent_values = {};
}

int AbstractValue::get_ID() {
  return(ID);
}

void AbstractValue::add_host_vector(RateVector* rv) {
  host_vectors.push_back(rv);
}

void AbstractValue::refresh_host_vectors() {
  for(auto it = host_vectors.begin(); it != host_vectors.end(); ++it) {
    (*it)->update_single_logLikelihood(ID);
  }
}

void AbstractValue::add_dependancy(AbstractValue* v) {
  dependent_values.push_back(v);
}

void AbstractValue::refresh() {
}

std::list<AbstractValue*> AbstractValue::get_dependancies() {
  return(dependent_values);
}

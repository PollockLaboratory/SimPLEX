#include <iostream>
#include <stdlib.h> //This gives rand.
#include <limits>

#include "CategoryFloat.h" 

// double inf = std::numeric_limits<double>::infinity();

CategoryFloat::CategoryFloat(std::string name, std::list<AbstractValue*> categories) : SampleableValue(name) {
  /*
   * The default constructor for the Continuous Float parameter class.
   */
  previous_value = 0.0;
}

// Utils
void CategoryFloat::printValue() {
  std::cout << "Category float - " << name << ": " << value << std::endl;
}

bool CategoryFloat::sample() {
  previous_value = value;
  fixedQ = false;

  value = 0.0;
  
  return(true);
}

double CategoryFloat::getValue() {
  return(value);
}

double CategoryFloat::getOldValue() {
  if(fixedQ) {
    std::cout << "Error: in CategoryFloat::getOldValue - already fixed." << std::endl;
    exit(EXIT_FAILURE);
  }
  return(previous_value);
}

void CategoryFloat::undo() {
  value = previous_value;
  previous_value = 0;
  fixedQ = true;
}

void CategoryFloat::fix() {
  previous_value = 0;
  fixedQ = true;
}


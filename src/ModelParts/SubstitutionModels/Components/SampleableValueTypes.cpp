#include "SampleableValueTypes.h"
#include <stdlib.h> //This gives rand.
#include <limits>

// CONTINUOUS FLOAT

double inf = std::numeric_limits<double>::infinity();

ContinuousFloat::ContinuousFloat(std::string name, int id, double initial_value = 0.0, double initial_std_dev = 1.0) : SampleableValue(name, id), value(initial_value), std_dev(initial_std_dev) {
  /*
   * The default constructor for the Continuous Float parameter class.
   */
  lower_bound = -inf;
  upper_bound = inf;

  previous_value = initial_value;
}

ContinuousFloat::ContinuousFloat(std::string name, int id, double initial_value = 0.0, double initial_std_dev = 1.0, double lower_bound = -inf) : SampleableValue(name, id), value(initial_value), std_dev(initial_std_dev), lower_bound(lower_bound) {
  /*
   * The default constructor for the Continuous Float parameter class - with bounds
   */
  upper_bound = inf;

  previous_value = initial_value;
}

ContinuousFloat::ContinuousFloat(std::string name, int id, double initial_value = 0.0, double initial_std_dev = 1.0, double lower_bound = -inf, double upper_bound = inf) : SampleableValue(name, id), value(initial_value), std_dev(initial_std_dev), lower_bound(lower_bound), upper_bound(upper_bound) {
  /*
   * The default constructor for the Continuous Float parameter class - with bounds
   */
  previous_value = initial_value;
}

void ContinuousFloat::print() {
  std::cout << "Continuous float - " << name << ": " << value << std::endl;
}

bool ContinuousFloat::sample() {
  previous_value = value;
  fixedQ = false;

  double r = ((rand() % 10000) / 10000.0) - 0.5;
  value = value + (r * std_dev);

  if(value < lower_bound) {
    value = 2*lower_bound - value;
  }

  if(value > upper_bound) {
    value = 2*upper_bound - value;
  }
  return(true);
}

const double& ContinuousFloat::getValue() {
  return(value);
}

const double& ContinuousFloat::getOldValue() {
  if(fixedQ) {
    std::cout << "Error: in ContinuousFloat::getOldValue - already fixed." << std::endl;
    exit(EXIT_FAILURE);
  }
  return(previous_value);
}

void ContinuousFloat::undo() {
  value = previous_value;
  previous_value = 0;
  fixedQ = true;
}

void ContinuousFloat::fix() {
  previous_value = 0;
  fixedQ = true;
}

void ContinuousFloat::refresh() {
}

std::ostream& operator<<(std::ostream& os, const ContinuousFloat& cf) {
  os << "[ContinuousFloat-" << cf.name << "]";
  return(os);
}

// CATEGORY FLOAT

DiscreteFloat::DiscreteFloat(std::string name, int id, RateCategories* categories) : SampleableValue(name, id) {
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

bool DiscreteFloat::sample() {
  prev_i = i;
  previous_value = value;
  fixedQ = false;

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
  
  return(true);
}

const double& DiscreteFloat::getValue() {
  return(value);
}

const double& DiscreteFloat::getOldValue() {
  if(fixedQ) {
    std::cout << "Error: in CategoryFloat::getOldValue - already fixed." << std::endl;
    exit(EXIT_FAILURE);
  }
  return(previous_value);
}

void DiscreteFloat::undo() {
  i = prev_i;
  value = previous_value;
  fixedQ = true;
}

void DiscreteFloat::fix() {
  fixedQ = true;
}

void DiscreteFloat::refresh(){
}

#include "ComponentTypes.h"

// RATE CATEGORIES

RateCategories::RateCategories(std::string name, float lower_b, float upper_b, int steps) : AbstractComponent(name), n(steps) {
  values = {};
  float step_size = (upper_b - lower_b) / steps;
  for(int i = 1; i <= steps; i++) {
    values.push_back(i*step_size);
  }
}

void RateCategories::refresh() {
}

void RateCategories::print() {
  std::cout << "Rate Categories - " << name << ": [ ";
  for(int i = 0; i < n; i++) {
    std::cout << values[i] << " ";
  }
  std::cout << "]" << std::endl;
}

float &RateCategories::operator[](int i) {
  return(values[i]);
}

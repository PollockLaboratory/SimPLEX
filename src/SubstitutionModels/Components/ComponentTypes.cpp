#include "ComponentTypes.h"

// RATE CATEGORIES

RateCategories::RateCategories(std::string name, int id, std::vector<float> categories) : AbstractComponent(name, id) {
  values = categories;
  n = values.size();
}

RateCategories::RateCategories(std::string name, int id, float lower_b, float upper_b, int steps) : AbstractComponent(name, id), n(steps) {
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

#ifndef ComponentTypes_h
#define ComponentTypes_h

#include "AbstractComponent.h"

#include <string>
#include <vector>
#include <list>
#include <iostream>

class RateCategories : public AbstractComponent {
 public:
  RateCategories(std::string name, std::vector<float> categories);
  RateCategories(std::string name, float lower_b, float upper_b, int steps);
  virtual void refresh();
  virtual void print();
  float &operator[](int);
  int n;
 private:
  std::vector<float> values;
};

#endif

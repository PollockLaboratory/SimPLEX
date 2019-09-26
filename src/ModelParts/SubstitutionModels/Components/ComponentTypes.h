#ifndef ComponentTypes_h
#define ComponentTypes_h

#include <string>
#include <vector>
#include <list>
#include <iostream>

#include "../../AbstractComponent.h"

class RateCategories : public AbstractComponent {
  // The vector representing the possible rate a discretely sampled value can be.
public:
  RateCategories(std::string name, int id, std::vector<float> categories);
  RateCategories(std::string name, int id, float lower_b, float upper_b, int steps);
  virtual void refresh();
  virtual void print();
  float &operator[](int);
  int n;
private:
  std::vector<float> values;
};

#endif

#ifndef FixedFloat_h_
#define FixedFloat_h_

#include "AbstractValue.h"
#include <string>

class FixedFloat : public AbstractValue {
 public:
  FixedFloat(std::string parameter_name, double);
  virtual double getValue();
  virtual double getOldValue();
  virtual void printValue();
  virtual void refresh();
 private:
  double value;
};

#endif

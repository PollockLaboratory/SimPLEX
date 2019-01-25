#ifndef ContinuousFloat_h_
#define ContinuousFloat_h_

#include <string>
#include <iostream>

#include "../AbstractValue.h"

class ContinuousFloat : public SampleableValue {
 public:
  ContinuousFloat(std::string, double, double);
  ContinuousFloat(std::string, double, double, double);
  ContinuousFloat(std::string, double, double, double, double);
  virtual void printValue();
  virtual bool sample();

  virtual double getValue();
  virtual double getOldValue();

  virtual void undo();
  virtual void fix();
  virtual void refresh();
 private:
  double value;
  double std_dev;

  double lower_bound;
  double upper_bound;

  double previous_value;	
};

#endif

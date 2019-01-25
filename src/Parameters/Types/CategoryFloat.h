#ifndef CategoryFloat_h_
#define CategoryFloat_h_

#include <string>
#include <list>
#include <iostream>

#include "../AbstractValue.h"

class CategoryFloat : public SampleableValue {
 public:
  CategoryFloat(std::string, std::list<AbstractValue*>);
  virtual void printValue();
  virtual bool sample();

  virtual double getValue();
  virtual double getOldValue();

  virtual void undo();
  virtual void fix();
 private:
  int i;
  double value;
  double previous_value;	
};

#endif

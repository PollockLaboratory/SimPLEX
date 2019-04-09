#ifndef SampleableValueTypes_h
#define SampleableValueTypes_h

#include "AbstractComponent.h"
#include "ComponentTypes.h"

#include <string>
#include <list>
#include <iostream>
#include <exception>

class ContinuousFloat : public SampleableValue {
 public:
  ContinuousFloat(std::string, int, double, double);
  ContinuousFloat(std::string, int, double, double, double);
  ContinuousFloat(std::string, int, double, double, double, double);

  virtual void print();
  virtual bool sample();

  virtual const double& getValue();
  virtual const double& getOldValue();

  virtual void undo();
  virtual void fix();
  virtual void refresh();

  friend std::ostream& operator<<(std::ostream&, const ContinuousFloat&);
 private:
  double value;
  double std_dev;

  double lower_bound;
  double upper_bound;

  double previous_value;	
};

class CategoryFloat : public SampleableValue {
 public:
  CategoryFloat(std::string, int id, RateCategories*);
  virtual void print();
  virtual bool sample();

  virtual const double& getValue();
  virtual const double& getOldValue();

  virtual void undo();
  virtual void fix();

  virtual void refresh();
 private:
  RateCategories* rc;
  int i;
  int prev_i;
  int n;
  double value;
  double previous_value;	
};

#endif

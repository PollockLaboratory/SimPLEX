#ifndef AbstractValueTypes_h
#define AbstractValueTypes_h

#include"AbstractComponent.h"
#include <string>
#include <list>
#include <iostream>
#include <exception>

class FixedFloat : public AbstractValue {
 public:
  FixedFloat(std::string parameter_name, double);
  virtual const double& getValue();
  virtual const double& getOldValue();
  virtual void print();
  virtual void refresh();
 private:
  double value;
};

class OutOfBoundsException: public std::exception {
 private:
  std::string m_error;
 public:
  OutOfBoundsException(std::string error_message) : m_error(error_message) {
  }
};

class VirtualSubstitutionRate : public AbstractValue {
 public:
  VirtualSubstitutionRate(std::string parameter_name, UniformizationConstant* unif);
  virtual const double& getValue();
  virtual const double& getOldValue();
  virtual void print();
  void refresh();
  void add_rate(AbstractValue*);
 private:
  UniformizationConstant* u;
  double value;
  double previous_value;
  std::list<AbstractValue*> dependent_rates;
};

#endif

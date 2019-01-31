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
  virtual double getValue();
  virtual double getOldValue();
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
  VirtualSubstitutionRate(std::string parameter_name, double unif);
  virtual double getValue();
  virtual double getOldValue();
  virtual void print();
  void refresh();
  void add_rate(AbstractValue*);
 private:
  double value;
  double previous_value;
  double u; //Uniformization constant.
  std::list<AbstractValue*> dependent_rates;
};

#endif

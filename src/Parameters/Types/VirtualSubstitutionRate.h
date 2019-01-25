#ifndef VirtualSubstitutionRate_h_
#define VirtualSubstitutionRate_h_

#include "AbstractValue.h"
#include <string>
#include <list>
#include <exception>

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
  virtual void printValue();
  void refresh();
  void add_rate(AbstractValue*);
 private:
  double value;
  double previous_value;
  double u; //Uniformization constant.
  std::list<AbstractValue*> dependent_rates;
};
#endif

#ifndef PARAMETERS_h
#define PARAMETERS_h

#include "../AbstractComponent.h"
//#include "../../IO/RawParameterTypes.h"

#include <string>
#include <list>
#include <iostream>
#include <exception>

// Sampleable

class ContinuousFloat : public SampleableValue {
public:
  ContinuousFloat(std::string, double, double);

  void print() override;
  sample_status sample() override;

  const double& get_value() override;
  const double& get_old_value() override;

  void undo() override;
  void fix() override;
  void refresh() override;

  std::string get_state_header() override;
  std::string get_state() override;

  std::string get_type() override;

  friend std::ostream& operator<<(std::ostream&, const ContinuousFloat&);
private:
  double value;
  double std_dev;

  double previous_value;	
};

class RateCategories : public AbstractComponent {
  // The vector representing the possible rate a discretely sampled value can be.
public:
  RateCategories(std::string name, std::vector<Valuable*> categories);
  void fix() override;
  void refresh() override;
  void print() override;

  std::string get_state_header() override;
  std::string get_state() override;

  std::string get_type() override;

  double operator[](int);
  int n;
private:
  std::vector<Valuable*> values;
};

class DiscreteFloat : public SampleableValue {
public:
  DiscreteFloat(std::string, RateCategories*);
  void print() override;
  sample_status sample() override;

  const double& get_value() override;
  const double& get_old_value() override;

  void undo() override;
  void fix() override;
  void refresh() override;

  std::string get_state_header() override;
  std::string get_state() override;

  std::string get_type() override;
private:
  RateCategories* rc;
  int i;
  int prev_i;
  int n;
  double value;
  double previous_value;	
};

// Dependent

class FixedFloat : public NonSampleableValue {
public:
  FixedFloat(std::string parameter_name, double);
  const double& get_value() override;
  const double& get_old_value() override;
  void print() override;
  std::string get_type() override;

  void fix() override;
  void refresh() override;

private:
  double value;
};

// Arithmatic

enum arith_op {
	       ADDITION,
	       SUBTRACTION,
	       MULTIPLICATION,
	       DIVISION
};

class Arithmatic : public NonSampleableValue{
public:
  static int idc;
  Arithmatic(arith_op, Valuable*, Valuable*);
  Arithmatic(std::string, arith_op, Valuable*, Valuable*);

  const double& get_value() override;
  const double& get_old_value() override;
  void print() override;

  void fix() override;
  void refresh() override;

  std::string get_type() override;
private:
  arith_op op;
  double (*arith_func)(double, double);
  double value;
  double previous_value;
  Valuable* v1;
  Valuable* v2;
  void setup(Valuable*, Valuable*);
};

// Dependency Groups

class DependencyGroup : public AbstractComponent {
public:
  DependencyGroup(std::string name);
  void set_parent(Valuable*);
  void print() override;

  void fix() override;
  void refresh() override;

  std::string get_state_header() override;
  std::string get_state() override;

  std::string get_type() override;
private:
  unsigned int n; // Number of parameters in group.
  Valuable* parent_parameter; // This is the parameter that all group elements reference when they are in the group.
};

class DependencyElement : public SampleableValue {
public:
  DependencyElement(std::string, SampleableValue*);
  void print() override;
  sample_status sample() override;

  const double& get_value() override;
  const double& get_old_value() override;

  void undo() override;
  void fix() override;
  void refresh() override;

  std::string get_state_header() override;
  std::string get_state() override;

  std::string get_type() override;
private:
  double value;
  DependencyGroup* group;
  SampleableValue* external_parameter; // The parameter that determines the value when not in a group.
  bool unified;
};

// Virtual Substitution rate.

class OutOfBoundsException: public std::exception {
private:
  std::string m_error;
public:
  OutOfBoundsException(std::string error_message) : m_error(error_message) {
  }
};

class VirtualSubstitutionRate : public NonSampleableValue {
public:
  VirtualSubstitutionRate(std::string name);
  const double& get_value() override;
  const double& get_old_value() override;
  void print() override;
  std::string get_type() override;

  void fix() override;
  void refresh() override;

  void set_u(Valuable* unif);
  void add_rate(Valuable* v);
private:
  Valuable* u;
  double value;
  double previous_value;
  std::list<Valuable*> dependent_rates;
};

#endif

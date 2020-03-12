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

  const double& getValue() override;
  const double& getOldValue() override;

  void undo() override;
  void fix() override;
  void refresh() override;

  double record_state(int gen, double l) override;
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

  double record_state(int gen, double l) override;
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

  const double& getValue() override;
  const double& getOldValue() override;

  void undo() override;
  void fix() override;
  void refresh() override;

  double record_state(int gen, double l) override;
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

class FixedFloat : public StaticValue {
 public:
  FixedFloat(std::string parameter_name, double);
  const double& getValue() override;
  const double& getOldValue() override;
  void print() override;

  void fix() override;
  void refresh() override;
  std::string get_type() override;
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

class VirtualSubstitutionRate : public StaticValue {
 public:
  VirtualSubstitutionRate(AbstractComponent* parameter, Valuable* unif);
  const double& getValue() override;
  const double& getOldValue() override;
  void print() override;
  std::string get_type() override;

  void fix() override;
  void refresh() override;
  void add_rate(Valuable* v);
 private:
  Valuable* u;
  double value;
  double previous_value;
  std::list<Valuable*> dependent_rates;
};

#endif

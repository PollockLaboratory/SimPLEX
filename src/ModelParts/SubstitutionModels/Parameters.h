#ifndef PARAMETERS_h
#define PARAMETERS_h

#include "../AbstractComponent.h"

#include <string>
#include <list>
#include <iostream>
#include <exception>

// Sampleable

class ContinuousFloat : public SampleableValue {
public:
  ContinuousFloat(std::string, double, double);
  ContinuousFloat(std::string, double, double, double);
  ContinuousFloat(std::string, double, double, double, double);

  virtual void print();
  virtual sample_status sample();

  virtual const double& getValue();
  virtual const double& getOldValue();

  virtual void undo();
  virtual void fix();
  virtual void refresh();

  virtual double record_state(int gen, double l);

  friend std::ostream& operator<<(std::ostream&, const ContinuousFloat&);
private:
  double value;
  double std_dev;

  double lower_bound;
  double upper_bound;

  double previous_value;	
};

class RateCategories : public AbstractComponent {
  // The vector representing the possible rate a discretely sampled value can be.
public:
  RateCategories(std::string name, std::vector<float> categories);
  RateCategories(std::string name, float lower_b, float upper_b, int steps);
  virtual void refresh();
  virtual void print();

  virtual double record_state(int gen, double l);

  float &operator[](int);
  int n;
private:
  std::vector<float> values;
};

class DiscreteFloat : public SampleableValue {
public:
  DiscreteFloat(std::string, RateCategories*);
  virtual void print();
  virtual sample_status sample();

  virtual const double& getValue();
  virtual const double& getOldValue();

  virtual void undo();
  virtual void fix();
  virtual void refresh();

  virtual double record_state(int gen, double l);
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

class VirtualSubstitutionRate : public StaticValue {
 public:
  VirtualSubstitutionRate(std::string parameter_name, Valuable* unif);
  virtual const double& getValue();
  virtual const double& getOldValue();
  virtual void print();
  void refresh();
  void add_rate(Valuable* v);
 private:
  Valuable* u;
  double value;
  double previous_value;
  std::list<Valuable*> dependent_rates;
};

#endif

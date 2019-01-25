#ifndef AbstractValue_h_
#define AbstractValue_h_

#include <string>
#include <vector>
#include <list>

// Given these current class definitions there cannot be a samplable parameter that is also dependent on other
// parameters - such as category parameters.

// Defined in RateVector.h
class RateVector;

class AbstractComponent {
 public:
  AbstractComponent(std::string name);

  void add_dependancy(AbstractComponent*);
  std::list<AbstractComponent*> get_dependancies();

  int get_ID();
  std::string get_name();

  virtual void refresh() = 0;
  virtual void printValue() = 0;
 protected:
  int ID;
  std::string name;
  std::list<AbstractComponent*> dependent_values; 
};

class AbstractValue : public AbstractComponent {
 public:
  AbstractValue(std::string name);
  virtual double getValue() = 0;
  virtual double getOldValue() = 0;
  virtual void printValue() = 0;
  // virtual void refresh() = 0;

  void add_host_vector(RateVector*);
  void refresh_host_vectors();

  std::list<RateVector*> host_vectors; // Pointers to the host RateVectors that a parameter sits within.
};

class SampleableValue : public AbstractValue {
 public:
  SampleableValue(std::string parameter_name) : AbstractValue(parameter_name) { fixedQ = true; }

  virtual bool sample() = 0; // If return true Metropolis Hasting Sample, else Gibbs.
  virtual void undo() = 0;
  virtual void fix() = 0;
  virtual void refresh() = 0;
 protected:
  bool fixedQ; //Indicates whether the state of this parameter is fixed or still being trialed.
};

#endif

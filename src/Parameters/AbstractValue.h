#ifndef AbstractValue_h_
#define AbstractValue_h_

#include <string>
#include <vector>
#include <list>

// Given these current class definitions there cannot be a samplable parameter that is also dependent on other
// parameters - such as category parameters.

// Defined in RateVector.h
class RateVector;

class AbstractValue {
 public:
  AbstractValue(std::string parameter_name);

  std::string name;
  int get_ID();
  virtual double getValue() = 0;
  virtual double getOldValue() = 0;
  virtual void printValue() = 0;

  virtual void refresh() = 0;

  void add_host_vector(RateVector*);
  void refresh_host_vectors();

  void add_dependancy(AbstractValue*);
  std::list<AbstractValue*> get_dependancies();
  std::list<RateVector*> host_vectors; // Pointers to the host RateVectors that a parameter sits within.
 protected:
  std::list<AbstractValue*> dependent_values;
 private:
  int ID;
};

class AbstractParameter : public AbstractValue {
 public:
  AbstractParameter(std::string parameter_name) : AbstractValue(parameter_name) { fixedQ = true; }

  virtual bool sample() = 0; // If return true Metropolis Hasting Sample, else Gibbs.
  virtual void undo() = 0;
  virtual void fix() = 0;
  void refresh() {}
 protected:
  bool fixedQ; //Indicates whether the state of this parameter is fixed or still being trialed.
};

#endif

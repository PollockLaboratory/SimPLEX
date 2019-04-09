#ifndef AbstractValue_h_
#define AbstractValue_h_

#include <string>
#include <vector>
#include <list>

// Given these current class definitions there cannot be a samplable parameter that is also dependent on other
// parameters - such as category parameters.

// Defined in RateVector.h
class RateVector;

struct rv_loc {
  // Rate Vector Locations
  RateVector* rv;
  int pos;
};

class AbstractComponent {
 public:
  AbstractComponent(std::string name, int id);

  void add_dependancy(AbstractComponent*);
  const std::list<AbstractComponent*>& get_dependancies();

  int get_ID();
  std::string get_name();

  virtual void refresh() = 0;
  virtual void print() = 0;
  std::list<rv_loc> host_vectors; // Abstract components do not use this but it is here to avoid dynamic_casting. Pointers to the host RateVectors that a parameter sits within.
  int ID;
  std::string name;
 protected:
  std::list<AbstractComponent*> dependent_values;
};

class AbstractValue : public AbstractComponent {
 public:
  AbstractValue(std::string name, int id);
  virtual const double& getValue() = 0;
  virtual const double& getOldValue() = 0;

  void add_host_vector(RateVector*, int);
  std::list<rv_loc> get_host_vectors();
};

class SampleableValue : public AbstractValue {
 public:
  SampleableValue(std::string parameter_name, int id) : AbstractValue(parameter_name, id) { fixedQ = true; }

  virtual bool sample() = 0; // If return true Metropolis Hasting Sample, else Gibbs.
  virtual void undo() = 0;
  virtual void fix() = 0;
  virtual void refresh() = 0;
 protected:
  bool fixedQ; //Indicates whether the state of this parameter is fixed or still being trialed.
};

class UniformizationConstant : public SampleableValue {
 public:
  UniformizationConstant();

  virtual void print();
  virtual bool sample();

  virtual const double& getValue();
  virtual const double& getOldValue();

  virtual void undo();
  virtual void fix();
  virtual void refresh();

  void add_VirtualSubstitutionRate(AbstractValue*);
  void set_initial();

  friend std::ostream& operator<<(std::ostream&, const UniformizationConstant&);
 private:
  double value;
  double previous_value;
  double threshold;
  double max_step;
  std::list<AbstractValue*> vsrs;
};

#endif

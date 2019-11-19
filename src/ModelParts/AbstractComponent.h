#ifndef AbstractValue_h_
#define AbstractValue_h_

#include <string>
#include <vector>
#include <list>

// Given these current class definitions there cannot be a samplable parameter that is also dependent on other
// parameters - such as category parameters.

// Defined in RateVector.h
class RateVector;

typedef struct sample_status {
  bool testp; // true if metropolis hasting sampling is required.
  bool updatedp; // true if the component has actually changed.
  bool full_recalculation; // True if full likelihood calculation is required.
} sample_status;

struct rv_loc {
  // Rate Vector Locations
  RateVector* rv;
  int pos;
};

class Valuable {
public:
  Valuable();
  virtual const double& getValue() = 0;
  virtual const double& getOldValue() = 0;
  void print();

  std::list<rv_loc> get_host_vectors();
  void add_host_vector(RateVector*, int);
  std::list<rv_loc> host_vectors; // Abstract components do not use this but it is here to avoid dynamic_casting. Pointers to the host RateVectors that a parameter sits within.
};

class AbstractComponent {
public:
  int ID;
  std::string name;
  AbstractComponent(std::string name);

  void add_dependancy(AbstractComponent*);
  const std::list<AbstractComponent*>& get_dependancies();

  int get_ID();
  std::string get_name();

  virtual void refresh() = 0;
  virtual void print() = 0;
  virtual double record_state(int gen, double l) = 0;
protected:
  std::list<AbstractComponent*> dependent_values;
};

class SampleableComponent : public AbstractComponent {
public:
  SampleableComponent(std::string name);
  virtual sample_status sample() = 0;
  virtual void undo() = 0;
  virtual void fix() = 0;
protected:
  bool fixedQ;
};

class SampleableValue : public SampleableComponent , public Valuable {
public:
  SampleableValue(std::string name);
};

class StaticValue : public AbstractComponent, public Valuable {
public:
  StaticValue(std::string name);
  virtual double record_state(int gen, double l);
};

class UniformizationConstant : public Valuable, public SampleableComponent {
public:
  UniformizationConstant(double);

  virtual void print();
  virtual sample_status sample();

  virtual const double& getValue();
  virtual const double& getOldValue();

  virtual void undo();
  virtual void fix();
  virtual void refresh();
  virtual double record_state(int gen, double l);

  void add_VirtualSubstitutionRate(Valuable*);
  void set_initial();

  friend std::ostream& operator<<(std::ostream&, const UniformizationConstant&);
private:
  double value;
  double previous_value;
  double threshold;
  double max_step;
  std::list<Valuable*> vsrs;
};

#endif

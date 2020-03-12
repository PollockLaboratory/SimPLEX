#ifndef AbstractValue_h_
#define AbstractValue_h_

#include <string>
#include <vector>
#include <list>
#include <set>

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
protected:
  // Dependencies -> this component -> dependents.
  std::list<AbstractComponent*> dependencies; // Components that this component depends on.
  std::list<AbstractComponent*> dependents; // Components that depend on this parameter.
  std::list<AbstractComponent*> refresh_list; // List of AbstractComponent that mush be refreshed when this AbstractComponent changes.
  std::list<Valuable*> valuable_dependents; // List of dependents of valuable type.

  std::list<AbstractComponent*> next_dependents(std::list<AbstractComponent*>, std::set<AbstractComponent*>&);
public:
  int ID;
  std::string name;
  AbstractComponent(std::string name);
  AbstractComponent(AbstractComponent* parameter);

  void add_dependancy(AbstractComponent*);
  const std::list<AbstractComponent*>& get_dependancies();

  void add_dependent(AbstractComponent*);
  const std::list<AbstractComponent*>& get_dependents();

  void setup_refresh_list();
  const std::list<AbstractComponent*>& get_refresh_list();
  const std::list<Valuable*>& get_valuable_dependents();

  int get_ID();
  std::string get_name();

  virtual void fix() = 0;
  virtual void refresh() = 0;
  virtual void print() = 0;
  virtual double record_state(int gen, double l) = 0;
  virtual std::string get_type() = 0;
};

class SampleableComponent : public AbstractComponent {
public:
  SampleableComponent(std::string name);
  virtual sample_status sample() = 0;
  virtual void undo() = 0;
protected:
  bool fixedQ;
};

// Constraints

class BaseConstraint {
public:
  virtual double get_value() const = 0;
  virtual std::string get_description() const = 0;
};

class FixedConstraint : public BaseConstraint {
private:
  double value;
public:
  FixedConstraint(double value);
  double get_value() const override;
  std::string get_description() const override;
};

class DynamicConstraint : public BaseConstraint {
private:
  Valuable* value;
public:
  DynamicConstraint(Valuable* value);
  double get_value() const override;
  std::string get_description() const override;
};

class SampleableValue : public SampleableComponent , public Valuable {
protected:
  BaseConstraint* lower_bound;
  BaseConstraint* upper_bound;
public:
  SampleableValue(std::string name);
  void set_lower_boundary(BaseConstraint*);
  void set_upper_boundary(BaseConstraint*);
};

class StaticValue : public AbstractComponent, public Valuable {
public:
  StaticValue(std::string name);
  StaticValue(AbstractComponent* parameter);
  double record_state(int gen, double l) override;
};

class UniformizationConstant : public Valuable, public SampleableComponent {
public:
  UniformizationConstant(double);

  void print() override;
  sample_status sample() override;

  const double& getValue() override;
  const double& getOldValue() override;

  void undo() override;
  void fix() override;
  void refresh() override;
  double record_state(int gen, double l) override;
  std::string get_type() override;

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

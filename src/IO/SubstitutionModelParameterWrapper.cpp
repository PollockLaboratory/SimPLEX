#include <limits>

#include "SubstitutionModelParameterWrapper.h"
#include "../ModelParts/SubstitutionModels/Parameters.h"
#include "LuaUtils.h"

double inf = std::numeric_limits<double>::infinity();
double neg_inf = -inf;

namespace IO {
  static Valuable* extract_float_or_parameter(const sol::object& val, std::string name, int i) {
    sol::optional<float> maybe_float = val.as<sol::optional<float>>();
    sol::optional<ParameterWrapper> maybe_parameter = val.as<sol::optional<ParameterWrapper>>();

    if(maybe_float) {
      return(new FixedFloat(name + "-" + std::to_string(i), maybe_float.value()));
    } else if(maybe_parameter) {
      Valuable* v = dynamic_cast<Valuable*>(maybe_parameter.value().parameter);
      if(v) {
	return(v);
      } else {
	std::cerr << "Error: parameter \'" << maybe_parameter.value().get_name() << "\' is not of type Valuable." << std::endl;
	  exit(EXIT_FAILURE);
      }
    } else {
      std::cerr << "Error: expecting elements of type float or parameter." << std::endl;
	exit(EXIT_FAILURE);
    }
  }

  // ParameterWrapper
  ParameterWrapper::ParameterWrapper(AbstractComponent* parameter) : name(parameter->get_name()) {
    this->parameter = parameter;
  }

  std::string ParameterWrapper::get_name() {
    return(this->name);
  }

  std::string ParameterWrapper::get_type() {
    // This isn't correct?
    return(parameter->get_type());
  }

  // Constraints.
  DynamicConstraint* create_bound(ParameterWrapper* param) {
    ContinuousFloat* cf = dynamic_cast<ContinuousFloat*>(param->parameter);
    if(cf == nullptr) {
      std::cerr << "Error: constraint parameter must be of type CONTINUOUS_FLOAT in boundary constraint." << std::endl;
      exit(EXIT_FAILURE);
    }

    return(new DynamicConstraint(cf));
  }
  
  void ParameterWrapper::set_lower_bound(ParameterWrapper* param) {
    if(this->get_type() != "CONTINUOUS_FLOAT") {
      std::cerr << "Error: setting the lower bound for parameter not of type CONTINUOUS_FLOAT." << std::endl;
      exit(EXIT_FAILURE);
    }

    DynamicConstraint* constraint = create_bound(param);

    SampleableValue* sv = dynamic_cast<SampleableValue*>(this->parameter);
    if(sv == nullptr) {
      std::cerr << "Error: unable to cast parameter to SampleableValue when setting lower bound." << std::endl;
      exit(EXIT_FAILURE); 
    }

    sv->set_lower_boundary(constraint);
  }

  void ParameterWrapper::set_upper_bound(ParameterWrapper* param) {
    if(this->get_type() != "CONTINUOUS_FLOAT") {
      std::cerr << "Error: setting the upper bound for parameter not of type CONTINUOUS_FLOAT." << std::endl;
      exit(EXIT_FAILURE);
    }

    DynamicConstraint* constraint = create_bound(param);

    SampleableValue* sv = dynamic_cast<SampleableValue*>(this->parameter);
    if(sv == nullptr) {
      std::cerr << "Error: unable to cast parameter to SampleableValue when setting lower bound." << std::endl;
      exit(EXIT_FAILURE); 
    }

    sv->set_upper_boundary(constraint);
  }


  // Creating parameters.
  ContinuousFloat* new_ContinuousFloat(std::string name, sol::table tbl) {
    double init = value_from_table<double>(tbl, "initial_value");
    double step_size = value_from_table<double>(tbl, "step_size");
    double lower_bound = tbl.get_or<double>("lower_bound", neg_inf);
    double upper_bound = tbl.get_or<double>("upper_bound", inf);

    ContinuousFloat* parameter = new ContinuousFloat(name, init, step_size);

    if(lower_bound != neg_inf) {
      parameter->set_lower_boundary(new FixedConstraint(lower_bound));
    }

    if(upper_bound != inf) {
      parameter->set_upper_boundary(new FixedConstraint(upper_bound));
    }

    return(parameter);
  }

  DiscreteFloat* new_DiscreteFloat(std::string name, sol::table tbl) {
    RateCategories* new_categories = dynamic_cast<RateCategories*>(value_from_table<ParameterWrapper>(tbl, "categories").parameter);
    //int initial_cat = tbl.get_or<int>("initial_category", 0);
    return(new DiscreteFloat(name, new_categories));
  }

  FixedFloat* new_FixedFloat(std::string name, sol::table tbl) {
    double value = value_from_table<double>(tbl, "value");
    return(new FixedFloat(name, value));
  }

  VirtualSubstitutionRate* new_VirtualSubstitutionRate(std::string name, sol::table) {
    return(new VirtualSubstitutionRate(name));
  }

  ParameterWrapper new_parameter(std::string name, std::string parameter_type, sol::table tbl) {
    AbstractComponent* param = nullptr;
    if(parameter_type == "continuous") {
      param = new_ContinuousFloat(name, tbl);
    } else if(parameter_type == "discrete") {
      param = new_DiscreteFloat(name, tbl);
    } else if(parameter_type == "fixed") {
      param = new_FixedFloat(name, tbl);
    } else if(parameter_type == "virtual") {
      param = new_VirtualSubstitutionRate(name, tbl);
    } else {
      std::cerr << "Error: " << parameter_type << " is not recognizes as a type of parameter." << std::endl;
      exit(EXIT_FAILURE);
    }
    return(ParameterWrapper(param));
  }

  ParameterWrapper new_categories(std::string name, sol::table tbl) {
    //std::vector<float> values = {};
    std::vector<Valuable*> values = {};
    for(auto kvp : tbl) {
      const sol::object& key = kvp.first;
      const sol::object& val = kvp.second;
      sol::optional<int> index = key.as<sol::optional<int>>();

      values.push_back(extract_float_or_parameter(val, name, index.value()));
    }
    return(new RateCategories(name, values));
  }

  // Arithmatic parameters.

  Valuable* extract_Valuable(ParameterWrapper param) {
    Valuable* val = dynamic_cast<Valuable*>(param.parameter);
    if(val == nullptr) {
      std::cerr << "Error: expecting a parameter in extract_Valuable" << std:: endl;
      exit(EXIT_FAILURE);
    } else {
      return(val);
    }
  }

  ParameterWrapper add_parameters(ParameterWrapper param1, ParameterWrapper param2) {
    Valuable* v1 = extract_Valuable(param1.parameter);
    Valuable* v2 = extract_Valuable(param2.parameter);
    AbstractComponent* param = new Arithmatic(ADDITION, v1, v2);
    return(ParameterWrapper(param));
  }

  ParameterWrapper named_add_parameters(std::string name, ParameterWrapper param1, ParameterWrapper param2) {
    Valuable* v1 = extract_Valuable(param1.parameter);
    Valuable* v2 = extract_Valuable(param2.parameter);
    AbstractComponent* param = new Arithmatic(name, ADDITION, v1, v2);
    return(ParameterWrapper(param));
  }

  ParameterWrapper subtract_parameters(ParameterWrapper param1, ParameterWrapper param2) {
    Valuable* v1 = extract_Valuable(param1.parameter);
    Valuable* v2 = extract_Valuable(param2.parameter);
    AbstractComponent* param = new Arithmatic(SUBTRACTION, v1, v2);
    return(ParameterWrapper(param));
  }

  ParameterWrapper named_subtract_parameters(std::string name, ParameterWrapper param1, ParameterWrapper param2) {
    Valuable* v1 = extract_Valuable(param1.parameter);
    Valuable* v2 = extract_Valuable(param2.parameter);
    AbstractComponent* param = new Arithmatic(name, SUBTRACTION, v1, v2);
    return(ParameterWrapper(param));
  }

  ParameterWrapper multiply_parameters(ParameterWrapper param1, ParameterWrapper param2) {
    Valuable* v1 = extract_Valuable(param1.parameter);
    Valuable* v2 = extract_Valuable(param2.parameter);
    AbstractComponent* param = new Arithmatic(MULTIPLICATION, v1, v2);
    return(ParameterWrapper(param));
  }

  ParameterWrapper named_multiply_parameters(std::string name, ParameterWrapper param1, ParameterWrapper param2) {
    Valuable* v1 = extract_Valuable(param1.parameter);
    Valuable* v2 = extract_Valuable(param2.parameter);
    AbstractComponent* param = new Arithmatic(name, MULTIPLICATION, v1, v2);
    return(ParameterWrapper(param));
  }

  ParameterWrapper divide_parameters(ParameterWrapper param1, ParameterWrapper param2) {
    Valuable* v1 = extract_Valuable(param1.parameter);
    Valuable* v2 = extract_Valuable(param2.parameter);
    AbstractComponent* param = new Arithmatic(DIVISION, v1, v2);
    return(ParameterWrapper(param));
  }

  ParameterWrapper named_divide_parameters(std::string name, ParameterWrapper param1, ParameterWrapper param2) {
    Valuable* v1 = extract_Valuable(param1.parameter);
    Valuable* v2 = extract_Valuable(param2.parameter);
    AbstractComponent* param = new Arithmatic(name, DIVISION, v1, v2);
    return(ParameterWrapper(param));
  }

  // DependencyGroups

  DependencyGroupWrapper::DependencyGroupWrapper(DependencyGroup* group) : name(group->get_name()) {
    this->group = group;
  }

  std::string DependencyGroupWrapper::get_name() {
    return(this->name);
  }

  DependencyGroupWrapper new_dependency_group(std::string name, sol::table) {
    DependencyGroup* group = new DependencyGroup(name);
    return(DependencyGroupWrapper(group));
  }
}

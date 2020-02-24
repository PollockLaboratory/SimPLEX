#include "RawParameterTypes.h"

#include "LuaUtils.h"
#include <limits>


namespace IO {
  double inf = std::numeric_limits<double>::infinity();
  double neg_inf = -inf;

  ParameterWrapper::ParameterWrapper(AbstractComponent* parameter) : name(parameter->name) {
    this->parameter = parameter;
  }

  std::string ParameterWrapper::get_name() {
    return(this->name);
  }

  std::string ParameterWrapper::get_type() {
    // This isn't correct.
    return(parameter->get_type());
  }

  ContinuousFloat* new_ContinuousFloat(std::string name, sol::table tbl) {
    double init = value_from_table<double>(tbl, "initial_value");
    double step_size = value_from_table<double>(tbl, "step_size");
    double lower_bound = tbl.get_or<double>("lower_bound", neg_inf);
    double upper_bound = tbl.get_or<double>("upper_bound", inf);
    return(new ContinuousFloat(name, init, step_size, lower_bound, upper_bound));
  }

  DiscreteFloat* new_DiscreteFloat(std::string name, sol::table tbl) {
    RateCategories* new_categories = dynamic_cast<RateCategories*>(value_from_table<ParameterWrapper>(tbl, "categories").parameter);
    return(new DiscreteFloat(name, new_categories));
  }

  ParameterWrapper new_parameter(std::string name, std::string parameter_type, sol::table tbl) {
    AbstractComponent* param = nullptr;
    if(parameter_type == "continuous") {
      param = new_ContinuousFloat(name, tbl);
    } else if(parameter_type == "discrete") {
      param = new_DiscreteFloat(name, tbl);
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
      sol::optional<float> maybe_float = val.as<sol::optional<float>>();
      sol::optional<ParameterWrapper> maybe_parameter = val.as<sol::optional<ParameterWrapper>>();
      if(maybe_float) {
	values.push_back(new FixedFloat(name + "-" + std::to_string(index.value()), maybe_float.value()));
	//values.push_back(maybe_float.value());
      } else if(maybe_parameter) {
	Valuable* v = dynamic_cast<Valuable*>(maybe_parameter.value().parameter);
	if(v) {
	  values.push_back(v);
	} else {
	  std::cerr << "Error: parameter \'" << maybe_parameter.value().get_name() << "'\ is not of type Valuable for RateCategories." << std::endl;
	  exit(EXIT_FAILURE);
	}
      } else {
	std::cerr << "Error: expecting elements of type float in categories \"" << name << "\"" << std::endl;
	exit(EXIT_FAILURE);
      }
    }
    return(new RateCategories(name, values));
  }
}

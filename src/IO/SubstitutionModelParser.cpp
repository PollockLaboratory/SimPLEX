#include "SubstitutionModelParser.h"

#include <string>
#include <iostream>

#include "../Environment.h"
#include "Files.h"
#include "LuaUtils.h"

extern IO::Files files;
extern Environment env;

sol::state lua;

namespace IO {

  // RAW PARAMETERS.
  raw_param::raw_param(std::string name, param_type t) : name(name), t(t) {
    static int i = 0;
    i++;
    ID = i;
  }

  std::ostream& operator<<(std::ostream& os, const raw_param& param) {
   os << "[" << param.ID << " " << param.name << "]";
   return(os);  
  }

  raw_param_wrapper::raw_param_wrapper(raw_param* ptr) : ptr(ptr) {
  }

  // Continuous Parameter
  raw_ContinuousFloat::raw_ContinuousFloat(std::string name, param_type t, sol::table tbl) : raw_param(name, t) {
    read_options_table(tbl);
  }

  void raw_ContinuousFloat::read_options_table(sol::table tbl) {
    init = value_from_table<double>(tbl, "initial_value");
    step_size = value_from_table<double>(tbl, "step_size");
  }

  // Discrete Parameter.
  raw_DiscreteFloat::raw_DiscreteFloat(std::string name, param_type t, sol::table tbl) : raw_param(name, t) {
    read_options_table(tbl);
  }

  void raw_DiscreteFloat::read_options_table(sol::table tbl) {
    sol::table cats_tbl = value_from_table<sol::table>(tbl, "categories");
    for(auto kvp : cats_tbl) {
      const sol::object& val = kvp.second;
      sol::optional<float> maybe_float = val.as<sol::optional<float>>();
      if(maybe_float) {
	categories.push_back(maybe_float.value());
      } else {
	std::cerr << "Error: expecting elements of type float in categories list in parameter \"" << name << "\"" << std::endl;
      }
    }
  }

  // RAW RATE VECTORS.

  raw_rate_vector::raw_rate_vector(std::string name, rv_use_class uc, std::list<raw_param*> rates) : name(name), uc(uc), rates(rates) {
    static int i = -1;
    i++;
    ID = i;
  }

  std::ostream& operator<<(std::ostream& os, const raw_rate_vector& rv) {
    os << rv.ID << "\t" << rv.name << " context: " << rv.uc.state << "\t[ ";
    for(auto it = rv.rates.begin(); it != rv.rates.end(); ++it) {
      os << *it << " ";
    }
    os << "]";
    return(os);
  }

  // RAW SUBSTITUTION MODEL.
  raw_substitution_model::raw_substitution_model() {
    states = {};
    ignore_states = {};
    params = {};
  }

  std::list<raw_param*> raw_substitution_model::get_parameters() {
	std::list<raw_param*> ret;
	for(auto it = rv_list.begin(); it != rv_list.end(); ++it) {
	  for(auto jt = (*it).rates.begin(); jt != (*it).rates.end(); ++jt) {
		ret.push_back(*jt);
	  }
	}
	return(ret);
  }

  const std::map<int, raw_param*>& raw_substitution_model::get_map_parameters() {
    return(params);
  }

  // Lua bindings
  void raw_substitution_model::set_states(sol::table tbl) {
    into_list(tbl, states);

    lua["states"]["count"] = states.size();

    int i = 1;
    for(auto it = states.begin(); it != states.end(); ++it) {
      lua["states"][i] = *it;
      i++;
    }
  }

  void raw_substitution_model::set_ignore_states(sol::table tbl) {
    into_list(tbl, ignore_states);
  }

  raw_param* raw_substitution_model::new_parameter(std::string name, std::string parameter_type, sol::table tbl) {
    raw_param* p;
    if(parameter_type == "continuous") {
      p = new raw_ContinuousFloat(name, FLOAT, tbl);
    } else if(parameter_type == "discrete") {
      p = new raw_DiscreteFloat(name, CATEGORY, tbl);
    } else {
      std::cerr << "Error: " << parameter_type << " is not recognizes as a type of parameter." << std::endl;
      exit(EXIT_FAILURE);
    }
    return(p);
  }

  std::list<int> get_positions(sol::table pos_tbl) {
    std::list<int> out = {};

    for(auto kvp : pos_tbl) {
      const sol::object& val = kvp.second;

      sol::optional<int> maybe_param = val.as<sol::optional<int>>();

      if(maybe_param) {
	int p = maybe_param.value();
	out.push_back(p);
      } else {
	std::cerr << "Error: expecting a integer in pos table for new rate vector."  << std::endl;
	exit(EXIT_FAILURE);
      }
    }
    return(out);
  }

  raw_rate_vector raw_substitution_model::new_rate_vector(std::string name, sol::table info_tbl, sol::table params_tbl) {
    std::string state = info_tbl["state"];
    std::list<int> possible_pos = get_positions(info_tbl["pos"]);
    rv_use_class uc = {state, possible_pos};
    std::list<raw_param*> rates = {};

    for(auto kvp : params_tbl) {
      const sol::object& val = kvp.second;

      sol::optional<raw_param_wrapper> maybe_param = val.as<sol::optional<raw_param_wrapper>>();

      if(maybe_param) {
	raw_param* p = maybe_param.value().ptr;
	rates.push_back(p);
      } else {
	std::cerr << "Error: expecting a Parameter."  << std::endl;
	exit(EXIT_FAILURE);
      }
    }

    raw_rate_vector rv(name, uc, rates);
    return(rv);
  }

  sol::table get_string_tbl(std::string option) {
    std::vector<std::string> vals = env.get_array<std::string>(option);
    sol::table tbl = lua.create_table();
    int i = 1;
    for(auto it = vals.begin(); it != vals.end(); ++it) {
      tbl[i] = *it;
      i++;
    }
    return(tbl);
  }

  void raw_substitution_model::add_rate_vector(raw_rate_vector rv) {
    rv_list.push_back(rv);
  }

  void raw_substitution_model::read_from_file(std::string file_name) {
    lua.open_libraries(sol::lib::base, sol::lib::table);

    //Main tables.
    auto model_table = lua["model"].get_or_create<sol::table>();
    model_table.set_function("set_name", [this](std::string name) -> void { this->name = name; });
    model_table.set_function("add_rate_vector", [this](raw_rate_vector rv) -> void { this->add_rate_vector(rv); });

    auto states_table = lua["states"].get_or_create<sol::table>();
    states_table.set_function("set", [this](sol::table tbl) -> void { this->set_states(tbl); });

    auto config_table = lua["config"].get_or_create<sol::table>();
    config_table.set_function("get_int", [](std::string key) -> int { return(env.get<int>(key)); });
    config_table.set_function("get_str", [](std::string key) -> std::string { return(env.get<std::string>(key)); });
    config_table.set_function("get_string_array", &get_string_tbl);

    lua.new_usertype<raw_param_wrapper>("Parameter",
					"new", [this](std::string name, std::string parameter_type, sol::table tbl) -> raw_param_wrapper {return(raw_param_wrapper(this->new_parameter(name, parameter_type, tbl)));});
    lua.new_usertype<raw_rate_vector>("RateVector",
    			      "new", [this](std::string name, sol::table info_tbl, sol::table param_tbl) -> raw_rate_vector { return(this->new_rate_vector(name, info_tbl, param_tbl)); });

    // Read the file.
    lua.script(files.read_all(file_name));

    // Tidy the output.
    for(auto it = rv_list.begin(); it != rv_list.end(); it++) {
      for(auto jt = (*it).rates.begin(); jt != (*it).rates.end(); ++jt) {
	params[(*jt)->ID] = *jt; 
      }
    } 
  }

  raw_substitution_model* read_substitution_model(std::string file_name) {
    raw_substitution_model* raw_model = new raw_substitution_model();
    raw_model->read_from_file(file_name);
    return(raw_model);
  }

  std::ostream& operator<<(std::ostream& os, const raw_substitution_model& sm) {
    os << "Raw Substitution Model: " << sm.name << std::endl;
    os << "States: ";
    for(auto it = sm.states.begin(); it != sm.states.end(); ++it) {
      os << "\"" << *it << "\" ";
    }
    os << std::endl;

    os << "Rate Vectors: " << std::endl;
    for(auto it = sm.rv_list.begin(); it != sm.rv_list.end(); ++it) {
      os << *it << std::endl;
    }
    return(os);
  }
}


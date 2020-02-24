#include <string>
#include <iostream>

#include "../Environment.h"
#include "Files.h"
#include "LuaUtils.h"

#include "SubstitutionModelParser.h"

extern IO::Files files;
extern Environment env;

sol::state lua;

namespace IO {
  // RAW RATE VECTORS.
  raw_rate_vector::raw_rate_vector(std::string name, rv_use_class uc, std::list<AbstractComponent*> new_rates) : name(name), uc(uc) {
    static int i = -1;
    i++;
    ID = i;
    this->rates = new_rates;
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
  }

  // Lua bindings & Util functions.
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

  raw_rate_vector new_rate_vector(std::string name, sol::table info_tbl, sol::table params_tbl) {
    std::string state = info_tbl["state"];
    std::list<int> possible_pos = get_positions(info_tbl["pos"]);
    rv_use_class uc = {state, possible_pos};
    std::list<AbstractComponent*> rates = {};

    for(auto kvp : params_tbl) {
      const sol::object& val = kvp.second;

      sol::optional<ParameterWrapper> maybe_param = val.as<sol::optional<ParameterWrapper>>();

      if(maybe_param) {
	rates.push_back(maybe_param.value().parameter);
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
    model_table.set_function("set_name", [this](std::string name) -> void {
					   this->name = name;
					 });
    model_table.set_function("add_rate_vector", [this](raw_rate_vector rv) -> void {
						  this->add_rate_vector(rv);
						});

    auto states_table = lua["states"].get_or_create<sol::table>();
    states_table.set_function("set", [this](sol::table tbl) -> void {
				       this->set_states(tbl);
				     });

    auto config_table = lua["config"].get_or_create<sol::table>();
    config_table.set_function("get_int", [](std::string key) -> int {
					   return(env.get<int>(key));
					 });

    config_table.set_function("get_str", [](std::string key) -> std::string {
					   return(env.get<std::string>(key));
					 });

    config_table.set_function("get_string_array", &get_string_tbl);

    lua.new_usertype<ParameterWrapper>("Parameter",
					"new", [](std::string name, std::string parameter_type, sol::table tbl) -> ParameterWrapper {
						 return(new_parameter(name, parameter_type, tbl));
					       },
					"name", &ParameterWrapper::get_name,
					"type", &ParameterWrapper::get_type);

    auto CatsTable = lua["Categories"].get_or_create<sol::table>();
    CatsTable.set_function("new", [](std::string name, sol::table tbl) -> ParameterWrapper {
				    return(new_categories(name, tbl));
    			  });
    
    lua.new_usertype<raw_rate_vector>("RateVector",
				      "new", [](std::string name, sol::table info_tbl, sol::table param_tbl) -> raw_rate_vector {
					       return(new_rate_vector(name, info_tbl, param_tbl));
					     });

    // Read the file.
    lua.script(files.read_all(file_name)); 
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


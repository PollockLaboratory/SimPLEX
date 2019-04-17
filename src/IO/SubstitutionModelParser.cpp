#include "SubstitutionModelParser.h"

#include <string>
#include <iostream>

#include "../Environment.h"
#include "Files.h"

extern IO::Files files;
extern Environment env;

sol::state lua;

namespace IO {

  raw_param::raw_param(std::string name, param_type t, double init) : name(name), t(t), init(init) {
    static int i = 0;
    i++;
    ID = i;
  }

  std::ostream& operator<<(std::ostream& os, const raw_param& param) {
    os << "[" << param.ID << " " << param.name << "]";
    return(os);  
  }

  raw_rate_vector::raw_rate_vector(std::string name, rv_use_class uc, std::list<raw_param> rates) : name(name), uc(uc), rates(rates) {
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
  }

  std::list<raw_param> raw_substitution_model::get_parameters() {
	std::list<raw_param> ret;
	for(auto it = rv_list.begin(); it != rv_list.end(); ++it) {
	  for(auto jt = (*it).rates.begin(); jt != (*it).rates.end(); ++jt) {
		ret.push_back(*jt);
	  }
	}
	return(ret);
  }

  // Lua bindings
  void raw_substitution_model::set_states(sol::table tbl) {
    //  std::list<std::string> states = {};
    for(auto kvp : tbl) {
      const sol::object& val = kvp.second;

      sol::optional<std::string> maybe_str = val.as<sol::optional<std::string>>();

      if(maybe_str) {
	std::string s = maybe_str.value();
	if(s == "-") {
	  std::cerr << "Error: the character \"-\" is reserved for gaps. Cannot be manually reassigned." << std::endl;
	  exit(EXIT_FAILURE);
	}
	states.push_back(s);
      } else {
	std::cerr << "Error: state is not String."  << std::endl;
	exit(EXIT_FAILURE);
      }
    }

    lua["states"]["count"] = states.size();

    int i = 1;
    for(auto it = states.begin(); it != states.end(); ++it) {
      lua["states"][i] = *it;
      i++;
    }
  }

  raw_param raw_substitution_model::new_parameter(std::string name, sol::table tbl) {
    std::string type = tbl["type"];
    sol::optional<double> initial_value = tbl["initial_value"];
    double init;
    if(initial_value) {
      init = initial_value.value();
    } else{
      std::cerr << "Error: initial_value is not set or is not a double." << std::endl;
      exit(EXIT_FAILURE);
    }

    param_type pt;
    if(type == "float") {
      pt = FLOAT;
    } else if(type == "int") {
      pt = INT;
    } else {
      std::cerr << "Error: " << type << " is not recognizes as a type of parameter." << std::endl;
      exit(EXIT_FAILURE);
    }
    raw_param p(name, pt, init);
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
    std::list<raw_param> rates = {};

    for(auto kvp : params_tbl) {
      const sol::object& val = kvp.second;

      sol::optional<raw_param> maybe_param = val.as<sol::optional<raw_param>>();

      if(maybe_param) {
	raw_param p = maybe_param.value();
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

  void raw_substitution_model::read_from_file(std::ifstream& lua_sm_in) {
    lua.open_libraries(sol::lib::base);

    //Main tables.
    auto model_table = lua["model"].get_or_create<sol::table>();
    model_table.set_function("set_name", [this](std::string name) -> void { this->name = name; });
    model_table.set_function("add_rate_vector", [this](raw_rate_vector rv) -> void { this->add_rate_vector(rv); });

    auto states_table = lua["states"].get_or_create<sol::table>();
    states_table.set_function("set", [this](sol::table tbl) -> void { this->set_states(tbl); });

    auto config_table = lua["config"].get_or_create<sol::table>();
    config_table.set_function("get_int", [](std::string key) -> int { return(env.get<int>(key)); });
    config_table.set_function("get_string_array", &get_string_tbl);

    lua.new_usertype<raw_param>("Parameter",
				"new", [this](std::string name, sol::table tbl) -> raw_param {return(this->new_parameter(name, tbl));});
    lua.new_usertype<raw_rate_vector>("RateVector",
				      "new", [this](std::string name, sol::table info_tbl, sol::table param_tbl) -> raw_rate_vector { return(this->new_rate_vector(name, info_tbl, param_tbl)); });

    // Read the file.
    std::stringstream buffer;
    buffer << lua_sm_in.rdbuf();

    lua.script(buffer.str());
  }

  raw_substitution_model* read_substitution_model(std::ifstream& lua_sm_in) {
    raw_substitution_model* raw_model = new raw_substitution_model();
    raw_model->read_from_file(lua_sm_in);
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


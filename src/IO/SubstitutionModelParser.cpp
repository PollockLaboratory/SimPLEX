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
    all_states["primary"] = {};
    ignore_states = {};
  }

  // STATES
  // Lua bindings & Util functions - with reading lua file.
  void raw_substitution_model::set_states(sol::table tbl) {
    all_states["primary"] = IO::readStates(extract_list(tbl));

    lua["States"]["count"] = all_states["primary"].size();

    int i = 1;
    for(auto it = all_states["primary"].begin(); it != all_states["primary"].end(); ++it) {
      lua["States"][i] = *it;
      i++;
    }
  }

  // Look into this - not sure the ignore states do anything.
  void raw_substitution_model::set_ignore_states(sol::table tbl) {
    ignore_states = extract_list(tbl);
  }

  void raw_substitution_model::add_hidden_state(std::string name, sol::table states, sol::table options) {
    std::list<std::string> new_states = {};
    try {
      new_states = IO::readStates(extract_list(states));
    } catch(IO::ParseException const &err) {
      throw IO::ParseException(std::string("error defining new hidden state: ") + err.what());
    }

    all_states[name] = new_states;

    // Output.
    if(options["sequences_output"] == nullptr) {
      states_seqs_output_files[name] = "";
    } else {
      states_seqs_output_files[name] = options["sequences_output"];
    }

    if(options["substitutions_output"] == nullptr) {
      states_subs_output_files[name] = "";
    } else {
      states_subs_output_files[name] = options["substitutions_output"];
    }
  }

  void raw_substitution_model::read_hidden_state_file(std::string hidden_state, std::string file_name) {
    auto it = all_states.find(hidden_state);
    if(it == all_states.end()) {
      std::cerr << "Error: attempting to read " << file_name << ", however \"" << hidden_state << "\" is not a recognized hidden state." << std::endl;
      exit(EXIT_FAILURE);
    }

    std::string fh = "hidden_state_"+hidden_state;
    files.add_file(fh, file_name, IOtype::INPUT);

    IO::RawAdvMSA MSA = {};
    try {
      MSA = readRawAdvMSA(files.read_all(fh), it->second);
    } catch(IO::ParseException const &err) {
      throw IO::ParseException(std::string("error reading ") + file_name + ": " + err.what());
    }

    hidden_states_data[hidden_state] = MSA;
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

  raw_rate_vector new_rate_vector(std::string name, sol::table info_tbl, sol::table params_tbl, std::map<std::string, std::list<std::string>> all_states) {

    // Find the state_set(domain) of the rate vector.
    // If not set, I will assume it is the primary domain.
    std::string state_domain;
    if(info_tbl["domain"] == nullptr) {
      state_domain = "primary";
    } else {
      state_domain = info_tbl["domain"];
      if(all_states.find(state_domain) == all_states.end()) {
	std::cerr << "Error: domain \"" << state_domain << "\" has not been specified." << std::endl;
	exit(EXIT_FAILURE);
      }
    }

    std::map<std::string, std::string> secondary_state = {};
    // Check Secondary States.
    for(auto it = all_states.begin(); it != all_states.end(); ++it) {
      if(it->first != state_domain) {
	if(info_tbl[it->first] == nullptr) {
	  std::cerr << "Error: specify for which state (in domain \"" << it->first << "\") rate vector " << name << " applies to." << std::endl;
	  exit(EXIT_FAILURE);
	}

	// Check valid state has been specified.
	secondary_state[it->first] = info_tbl[it->first];
      }
    }

    // Check the state is specified in the domain.
    std::set<std::string> state_set(all_states[state_domain].begin(), all_states[state_domain].end());
    std::string state = info_tbl["state"];
    if(state_set.find(state) == state_set.end()) {
      std::cerr << "Error: " << state << " is not specified within the " << state_domain << " domain." << std::endl;
      exit(EXIT_FAILURE);		 
    }

    std::list<int> possible_pos = get_positions(info_tbl["pos"]);
    rv_use_class uc = {state_domain, state, possible_pos, secondary_state};
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

  sol::table get_string_tbl(std::string option) { // from the Environment.
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
    lua.open_libraries(sol::lib::base, sol::lib::table, sol::lib::string);

    // Main Lua tables.
    // MODEL
    auto model_table = lua["Model"].get_or_create<sol::table>();
    model_table.set_function("set_name", [this](std::string name) -> void {
					   this->name = name;
					 });

    model_table.set_function("add_rate_vector", [this](raw_rate_vector rv) -> void {
						  this->add_rate_vector(rv);
						});

    // STATES
    auto states_table = lua["States"].get_or_create<sol::table>();
    states_table.set_function("set", [this](sol::table tbl) -> void {
				       this->set_states(tbl);
				     });

    states_table.set_function("new_hidden", [this](std::string name, sol::table states, sol::table options) -> void {
					      this->add_hidden_state(name, states, options);
					    });

    // DATA - tools to incorperate additional data into the state, primarily hidden states.
    auto data_table = lua["Data"].get_or_create<sol::table>();

    data_table.set_function("load_hidden_state", [this](std::string hidden_state, std::string file_name) -> void {
						   read_hidden_state_file(hidden_state, file_name);
						   
						 });

        // CONFIGURATION
    auto config_table = lua["Config"].get_or_create<sol::table>();
    config_table.set_function("get_int", [](std::string key) -> int {
					   return(env.get<int>(key));
					 });

    config_table.set_function("get_float", [](std::string key) -> double {
					     return(env.get<double>(key));
					   });

    config_table.set_function("get_str", [](std::string key) -> std::string {
					   return(env.get<std::string>(key));
					 });

    config_table.set_function("get_string_array", &get_string_tbl);

    // PARAMETERS
    lua.new_usertype<ParameterWrapper>("Parameter",
				       "new", [](std::string name, std::string parameter_type, sol::table tbl) -> ParameterWrapper {
						 return(new_parameter(name, parameter_type, tbl));
					       },
				       "name", &ParameterWrapper::get_name,
				       "type", &ParameterWrapper::get_type,
				       "set_lower_bound", &ParameterWrapper::set_lower_bound,
				       "set_upper_bound", &ParameterWrapper::set_upper_bound,
				       "add", add_parameters,
				       "named_add", named_add_parameters,
				       "subtract", subtract_parameters,
				       "named_subtract", named_subtract_parameters,
				       "multiply", multiply_parameters,
				       "named_multiply", named_multiply_parameters,
				       "divide", divide_parameters,
				       "named_divide", named_divide_parameters);

    // Not fully implimented at all.
    lua.new_usertype<DependencyGroupWrapper>("DependencyGroup",
					     "new", [](std::string name, sol::table tbl) -> DependencyGroupWrapper {
						      return(new_dependency_group(name, tbl));
						    },
					     "name", &DependencyGroupWrapper::get_name);

    // I think this can be tidied up a bit.
    auto CatsTable = lua["Categories"].get_or_create<sol::table>();
    CatsTable.set_function("new", [](std::string name, sol::table tbl) -> ParameterWrapper {
				    return(new_categories(name, tbl));
    			  });
    
    lua.new_usertype<raw_rate_vector>("RateVector",
				      "new", [this](std::string name, sol::table info_tbl, sol::table param_tbl) -> raw_rate_vector {
					       return(new_rate_vector(name, info_tbl, param_tbl, all_states));
					     });

    // Read the file.
    try {
      lua.script(files.read_all(file_name));
    } catch(const sol::error &err) {
      std::cerr << "Error when reading lua script specifying substitution model:" << std::endl;
      std::cerr << err.what() << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  const std::list<std::string> raw_substitution_model::get_states() {
    return(all_states["primary"]);
  }

  const std::map<std::string, std::list<std::string>> raw_substitution_model::get_all_states() {
    return(all_states);
  }

  const std::list<std::string> raw_substitution_model::get_ignore_states() {
    return(ignore_states);
  }

  const std::list<raw_rate_vector> raw_substitution_model::get_rate_vector_list() {
    return(rv_list);
  }

  const IO::RawAdvMSA raw_substitution_model::get_hidden_states_data(std::string name) {
    return(hidden_states_data[name]);
  }

  raw_substitution_model* read_substitution_model(std::string file_name) {
    raw_substitution_model* raw_model = new raw_substitution_model();
    raw_model->read_from_file(file_name);
    return(raw_model);
  }

  std::ostream& operator<<(std::ostream& os, raw_substitution_model& sm) {
    os << "Raw Substitution Model: " << sm.name << std::endl;
    os << "States: ";
    for(auto it = sm.all_states["primary"].begin(); it != sm.all_states["primary"].end(); ++it) {
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


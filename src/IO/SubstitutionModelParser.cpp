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
  raw_rate_vector::raw_rate_vector(std::string name, RVScope uc, std::list<AbstractComponent*> new_rates) : name(name), uc(uc) {
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
  }

  // STATES
  void raw_substitution_model::add_state(std::string name, sol::table states, sol::table options) {
    std::list<std::string> new_states = {};
    try {
      new_states = IO::readStates(extract_list(states));
    } catch(IO::ParseException const &err) {
      throw IO::ParseException(std::string("error defining new hidden state: ") + err.what());
    }

    all_states[name] = new_states;
    lua["States"][name] = new_states;

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

  std::optional<std::list<std::string>> get_states_list(std::map<std::string, std::list<std::string>> all_states, std::string state_domain) {
    auto it = all_states.find(state_domain);
    if(it == all_states.end()) return std::nullopt;
    return std::make_optional<std::list<std::string>>(it->second);
  }

  const IO::RawMSA* read_input_MSA(std::string state_domain, std::string file_name, std::list<std::string> states) {
    std::string fh = "state_"+state_domain;
    files.add_file(fh, file_name, IOtype::INPUT);

    IO::RawMSA *MSA = new IO::RawMSA();
    try {
      *MSA = readRawMSA(files.read_all(fh), states);
    } catch(IO::ParseException const &err) {
      throw IO::ParseException(std::string("error reading ") + file_name + ": " + err.what());
    }

    return MSA;
  }
  
  void raw_substitution_model::load_dynamic_from_file(std::string state_domain, std::string file_name) {
    auto opt_states = get_states_list(this->all_states, state_domain);
    if(not opt_states.has_value()) {
      std::cerr << "Error: attempting to read " << file_name << ", however \"" << state_domain << "\" is not a recognized hidden state." << std::endl;
      exit(EXIT_FAILURE);
    }

    StateData sd = { StateData::Tag::DYNAMIC, {} };
    sd.data.dynamic = read_input_MSA(state_domain, file_name, opt_states.value());

    this->state_data[state_domain] = sd;
  }

  void raw_substitution_model::load_site_static_from_file(std::string state_domain, std::string file_name) {
    auto opt_states = get_states_list(this->all_states, state_domain);
    if(not opt_states.has_value()) {
      std::cerr << "Error: attempting to read " << file_name << ", however \"" << state_domain << "\" is not a recognized hidden state." << std::endl;
      exit(EXIT_FAILURE);
    }

    StateData sd = { StateData::Tag::SITE_STATIC, {} };
    sd.data.site_static = read_input_MSA(state_domain, file_name, opt_states.value());

    this->state_data[state_domain] = sd;
  }

  void raw_substitution_model::generate_uniform_data(std::string domain) {
    std::cout << "Implement me" << std::endl;
    exit(EXIT_FAILURE);
    // Get the states.
    //std::list<std::string> states = all_states[domain];

    // Find a MSA that is filled with external data.
    // This template msa is used to find the location of gaps.
    //IO::RawMSA template_msa;
    //for(auto it = state_data.begin(); it != state_data.end(); ++it) {
    //  std::cout << it->first << " " << it->second->n << " " << it->second->cols << std::endl;
    //  if(it->first == domain) {
    //    std::cerr << "Error: cannot create a uniform prior for " << domain << ", as data for that domain already exists." << std::endl;
    //    exit(EXIT_FAILURE);
    //  } else {
    //    template_msa = it->second;
    //    break;
    //  }
    //}

    //IO::RawMSA msa = IO::createUniformPrior(states, template_msa);

    //IO::printRawAdvMSA(msa);

    //state_data[domain] = msa;
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
    // If not set, will be assumed its name is "primary".
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
      std::cerr << "Error: the use cases for the rate vector " << name << "are not valid." << std::endl;
      std::cerr << "The state \'" << state << "\' is not specified within the " << state_domain << " domain." << std::endl;
      exit(EXIT_FAILURE);		 
    }

    RVScope uc = {state_domain, state, secondary_state};
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
    lua.open_libraries(sol::lib::base, sol::lib::table, sol::lib::string, sol::lib::io, sol::lib::math, sol::lib::os);

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

    states_table.set_function("new", [this](std::string name, sol::table states, sol::table options) -> void {
					      this->add_state(name, states, options);
					    });

    // DATA - tools to incorperate additional data into the state, primarily hidden states.
    auto data_table = lua["Data"].get_or_create<sol::table>();

    data_table.set_function("load_state",
                            [this](std::string state, std::string file_name) -> void {
                              load_dynamic_from_file(state, file_name); 
                            });

    data_table.set_function("load_site_static_state",
                            [this](std::string state, std::string file_name) -> void {
                              load_site_static_from_file(state, file_name); 
                            }); 

    data_table.set_function("generate_uniform_data",
                            [this](std::string state) -> void {
                              std::cout << "Creating a uniform prior: " << state << std::endl;
                              generate_uniform_data(state);
                            });

    // CONFIGURATION
    auto config_table = lua["Config"].get_or_create<sol::table>();
    config_table.set_function("get_bool", [](std::string key) -> int {
					    return(env.get<bool>(key));
					  });
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

    config_table.set_function("get_root_directory", []() -> std::string {
						      return(files.get_root_directory());
						    });

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
    //lua.new_usertype<DependencyGroupWrapper>("DependencyGroup",
		//			     "new", [](std::string name, sol::table tbl) -> DependencyGroupWrapper {
		//				      return(new_dependency_group(name, tbl));
		//				    },
		//			     "name", &DependencyGroupWrapper::get_name);

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

  const std::map<std::string, std::list<std::string>> raw_substitution_model::get_all_states() {
    return(all_states);
  }

  const std::list<raw_rate_vector> raw_substitution_model::get_rate_vector_list() {
    return(rv_list);
  }

  const StateData raw_substitution_model::get_state_data(std::string name) {
    return(state_data[name]);
  }

  raw_substitution_model* read_substitution_model(std::string file_name) {
    raw_substitution_model* raw_model = new raw_substitution_model();
    raw_model->read_from_file(file_name);
    return(raw_model);
  }
}


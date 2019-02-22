#include "CustomModel.h"

#include <iostream>
#include <cstdlib>
#include <list>

#include "../../Environment.h"
#include "../../IO.h"

#include "sol2/sol.hpp"

extern Environment env;
extern IO::Files files;

CustomModel::CustomModel() : SubstitutionModel() {
}

void CustomModel::set_name(std::string n) {
  name = n;
}

void CustomModel::set_states(sol::table tbl) {
  std::list<std::string> states = {};
  for(auto kvp : tbl) {
    const sol::object& key = kvp.first;
    const sol::object& val = kvp.second;

    sol::optional<std::string> maybe_str = val.as<sol::optional<std::string>>();

    if(maybe_str) {
      std::string s = maybe_str.value();
      if(s == "-") {
	std::cerr << "Error: the character \"-\" is reserved for gaps. Cannot be manually reassigned." << std::endl;
	exit(EXIT_FAILURE);
      }
      add_state(s);
    } else {
      std::cerr << "Error: state is not String."  << std::endl;
      exit(EXIT_FAILURE);
    }
  }
}

RateVector* lua_RateVector_cstr(std::string name, int state, sol::table tbl) {

  std::vector<AbstractValue*> rates = {};
  for(auto kvp : tbl) {
    const sol::object& key = kvp.first;
    const sol::object& val = kvp.second;

    sol::optional<AbstractValue*> maybe_val = val.as<sol::optional<AbstractValue*>>();

    if(maybe_val) {
      AbstractValue* v = maybe_val.value();
      rates.push_back(v);
    } else {
      std::cerr << "Error: value is not AbstractValue."  << std::endl;
      exit(EXIT_FAILURE);
    }
  }
  // State decreased by one to reflect differance in indexing between lua and C++.
  return(new RateVector(name, state - 1, rates));
}

CategoryFloat* lua_CategoryFloat_cstr(std::string name, sol::table tbl) {
  std::vector<float> categories = {};

  for(auto kvp : tbl) {
    const sol::object& key = kvp.first;
    const sol::object& val = kvp.second;

    sol::optional<float> maybe_val = val.as<sol::optional<float>>();

    if(maybe_val) {
      float v = maybe_val.value();
      categories.push_back(v);
    } else {
      std::cerr << "Error in constructing categories set: value is not float."  << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  RateCategories* rc = new RateCategories(name, categories);
  return(new CategoryFloat(name, rc));
}

void CustomModel::Initialize(int number_of_sites, std::vector<std::string> states) {
}

void CustomModel::Initialize() {
  float u = env.u;

  sol::state lua;
  lua.open_libraries(sol::lib::base);

  // Main tables.
  auto model_table = lua["model"].get_or_create<sol::table>();
  model_table.set_function("set_name", [this](std::string name) -> void { this->set_name(name); });
  model_table.set_function("add_rate_vector", [this](RateVector* rv) -> void { this->add_rate_vector(rv); });

  auto states_table = lua["states"].get_or_create<sol::table>();
  states_table.set_function("set", [this](sol::table tbl) -> void { this->set_states(tbl); });
  states_table.set_function("print", [this]() -> void { this->print_states(); });
  states_table.set_function("to_int", [this](std::string s) -> int { return(this->states.state_to_int[s] + 1); });
  states_table.set_function("to_str", [this](int i) -> std::string { return(this->states.int_to_state[i - 1]); });

  // Defining Usertypes - differant parameter types.
  lua.new_usertype<VirtualSubstitutionRate>("VirtualSubstitutionRate",
					    sol::base_classes, sol::bases<AbstractValue>(),
					    "new", [](std::string name) -> VirtualSubstitutionRate* {
					      return(new VirtualSubstitutionRate(name, env.u));},
					    "add_rate", &VirtualSubstitutionRate::add_rate);

  lua.new_usertype<ContinuousFloat>("ContinuousFloat",
				    sol::base_classes, sol::bases<AbstractComponent, AbstractValue, SampleableValue>(),
				    "new", sol::overload([](std::string name, double i, double j) -> ContinuousFloat* { return(new ContinuousFloat(name, i, j)); },
							 [](std::string name, double i, double j, double k) -> ContinuousFloat* { return(new ContinuousFloat(name, i, j, k)); },
							 [](std::string name, double i, double j, double k, double l) -> ContinuousFloat* { return(new ContinuousFloat(name, i, j, k, l)); }));

  lua.new_usertype<CategoryFloat>("CategoryFloat",
				  sol::base_classes, sol::bases<AbstractComponent, AbstractValue, SampleableValue>(),
				  "new", &lua_CategoryFloat_cstr);
				  

  lua.new_usertype<FixedFloat>("FixedFloat",
			       sol::base_classes, sol::bases<AbstractComponent, AbstractValue>(),
			       "new", [](std::string name, double value) -> FixedFloat* { return(new FixedFloat(name, value)); });

  lua.new_usertype<RateVector>("RateVector",
			       "new", &lua_RateVector_cstr);

  // Read the Lua file.
  files.add_file("lua_model", env.get("custom_model"), IOtype::INPUT);
  std::cout << "Custom model.\n== Reading from file: " << env.get("custom_model") << " ==" << std::endl;

  std::ifstream lua_file = files.get_ifstream("lua_model");

  std::stringstream buffer;
  buffer << lua_file.rdbuf();

  lua.script(buffer.str());

  std::cout << "== Finished reading Custom model ==" << std::endl;

  finalize();
}

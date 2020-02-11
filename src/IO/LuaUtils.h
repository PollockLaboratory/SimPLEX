#ifndef LuaUtils_h_
#define LuaUtils_h_

#include <string>
#include <iterator>

template <class T>
static T value_from_table(sol::table tbl, std::string name) {
  sol::optional<T> opt_val = tbl[name];
  if(not opt_val) {
    std::cerr << "Error: " << name << "has not been correctly specified." << std::endl;
    exit(EXIT_FAILURE);
  } else {
    return(opt_val.value());
  }
}

void into_list(sol::table tbl, std::list<std::string>& list) {
  for(auto kvp : tbl) {
      const sol::object& val = kvp.second;

      sol::optional<std::string> maybe_str = val.as<sol::optional<std::string>>();

      if(maybe_str) {
	std::string s = maybe_str.value();
	if(s == "-") {
	  std::cerr << "Error: the character \"-\" is reserved for gaps. Cannot be manually reassigned." << std::endl;
	  exit(EXIT_FAILURE);
	}
	list.push_back(s);
      } else {
	std::cerr << "Error: state is not String."  << std::endl;
	exit(EXIT_FAILURE);
      }
    }
}

#endif

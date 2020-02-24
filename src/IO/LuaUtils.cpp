#include "LuaUtils.h"

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

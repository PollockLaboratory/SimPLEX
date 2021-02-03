#include "LuaUtils.h"

std::list<std::string> extract_list(sol::table tbl) {
  std::list<std::string> ret = {};
  for(auto kvp : tbl) {
    const sol::object& val = kvp.second;
    sol::optional<std::string> maybe_str = val.as<sol::optional<std::string>>();

    if(maybe_str) {
      std::string s = maybe_str.value();
      ret.push_back(s);
    } else {
      std::cerr << "Error: expecting list of strings." << std::endl;
      exit(EXIT_FAILURE);
    }
  }
  return(ret);
}

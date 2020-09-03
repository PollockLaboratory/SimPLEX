#ifndef LuaUtils_h_
#define LuaUtils_h_

#include <string>
#include <list>
#include <iterator>
#include <iostream>

#include "sol2/sol.hpp"

template <class T>
static T value_from_table(sol::table tbl, std::string name) {
  sol::optional<T> opt_val = tbl[name];
  if(not opt_val) {
    std::cerr << "Error: \'" << name << "\' has not been correctly set in Lua table." << std::endl;
    exit(EXIT_FAILURE);
  } else {
    return(opt_val.value());
  }
}

//template <class T>
//static T possible_value_from_table(sol::table tbl, std::string name) {
//  T ret = tbl.get(name, NULL);
//  return(ret);
//}

void into_list(sol::table tbl, std::list<std::string>& list);

#endif

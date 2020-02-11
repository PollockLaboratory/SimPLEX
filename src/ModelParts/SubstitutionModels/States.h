#ifndef States_h_
#define States_h_

#include <map>
#include <string>
#include <set>

struct States {
  int n; // Careful about indels.
  std::set<std::string> possible;
  std::map<std::string, int> state_to_int;
  std::map<int, std::string> int_to_state;
};

#endif

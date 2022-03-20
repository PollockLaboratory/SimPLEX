#ifndef States_h_
#define States_h_

#include <map>
#include <string>
#include <set>
#include <list>

typedef signed char state_element;

struct States {
  unsigned int n; // Careful about indels.
  std::set<std::string> possible;
  std::map<std::string, state_element> state_to_int;
  std::map<state_element, std::string> int_to_state;
};

void print_States(States);

States add_to_States(States, std::string);
States create_States(std::list<std::string>);

#endif

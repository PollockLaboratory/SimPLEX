#include "States.h"

#include <iostream>

void print_States(States states) {
  std::cout << "States: [ ";
 for(auto it = states.state_to_int.begin(); it != states.state_to_int.end(); ++it) {
   std::cout << it->first << ":" << it->second << " ";
 }
 std::cout << "]- n = " << states.possible.size() << std::endl;
}

States add_to_States(States states, std::string s) {
  states.possible.insert(s);
  states.state_to_int[s] = states.n;
  states.int_to_state[states.n] = s;
  states.n++;

  return(states);
}

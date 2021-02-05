#include "States.h"

#include <iostream>

void print_States(States states) {
  std::cout << "[ ";
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

States create_States(std::set<std::string> input_states) {
  std::map<std::string, int> state_to_int = {};
  std::map<int, std::string> int_to_state = {};

  int i = 0;
  for(auto it = input_states.begin(); it != input_states.end(); ++it) {
    state_to_int[*it] = i;
    int_to_state[i] = *it;
    i++;
  }

  States states = {input_states.size(), input_states, state_to_int, int_to_state};
  return(states);
}

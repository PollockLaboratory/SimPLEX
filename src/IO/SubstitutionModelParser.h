#ifndef SubstitutionModelParser_h_
#define SubstitutionModelParser_h_

#define SOL_EXCEPTIONS_SAFE_PROPAGATION

#include <string>
#include <iostream>
#include <list>
#include <map>

#include "sol2/sol.hpp"
#include "SubstitutionModelParameterWrapper.h"
#include "SequencesParser.h"

namespace IO {
  // RATE VECTOR
  struct rv_use_class {
    // Structure describing where a rate vector can apply.
    std::string state;
    std::list<int> pos;
  };

  class raw_rate_vector {
  public:
    raw_rate_vector(std::string, rv_use_class, std::list<AbstractComponent*>);
    int ID;
    std::string name;
    rv_use_class uc;
    std::list<AbstractComponent*> rates;
    friend std::ostream& operator<<(std::ostream&, const IO::raw_rate_vector&);
  };
  
  class raw_substitution_model {
  public:
    std::string name;
    raw_substitution_model();
    void read_from_file(std::string file_name);
    friend std::ostream& operator<<(std::ostream&, const IO::raw_substitution_model&);

    const std::set<std::string> get_states();
    const std::list<std::string> get_ignore_states();
    const std::list<raw_rate_vector> get_rate_vector_list();

    const std::map<std::string, std::set<std::string>> get_hidden_states();
  private:
    std::list<std::string> ignore_states;
    std::set<std::string> states;
    std::list<raw_rate_vector> rv_list;

    void set_states(sol::table tbl);
    void set_ignore_states(sol::table tbl);
    void add_rate_vector(raw_rate_vector rv);

    // Hidden states.
    std::map<std::string, std::set<std::string>> hidden_states;
    std::map<std::string, IO::RawAdvMSA> hidden_states_data;

    void add_hidden_state(std::string name, sol::table states, sol::table options);
    void read_hidden_state_file(std::string hidden_state, std::string file_name);
  };

  raw_substitution_model* read_substitution_model(std::string file_name);
}

#endif

#ifndef SubstitutionModelParser_h_
#define SubstitutionModelParser_h_

#include <string>
#include <iostream>
#include <list>
#include <map>

#include "sol2/sol.hpp"
#include "RawParameterTypes.h"

namespace IO {
  // RATE VECTOR
  struct rv_use_class {
    // Structure describing where a rate vector can apply.
    std::string state;
    std::list<int> pos;
  };

  class raw_rate_vector {
  public:
    //raw_rate_vector(std::string, rv_use_class, std::list<raw_Parameter*>);
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
    std::list<std::string> states;
    std::list<std::string> ignore_states;
    std::list<raw_rate_vector> rv_list;
    //std::list<AbstractComponent*> get_parameters();
    //const std::map<int, raw_Parameter*>& get_map_parameters();
  private:
    //void fill_parameter_list(raw_Parameter*);
    //std::map<int, raw_Parameter*> params;
    void set_states(sol::table tbl);
    void set_ignore_states(sol::table tbl);
    void add_rate_vector(raw_rate_vector rv);
  };

  raw_substitution_model* read_substitution_model(std::string file_name);
}

#endif

#ifndef SubstitutionModelParser_h_
#define SubstitutionModelParser_h_

#include <string>
#include <iostream>
#include <list>

#include "sol2/sol.hpp"

namespace IO {
  enum param_type { INT, FLOAT };

  class raw_param {
  public:
    raw_param(std::string, param_type, double);
    int ID;
    std::string name;
    param_type t;
    double init;
    friend std::ostream& operator<<(std::ostream&, const IO::raw_param&);
  };

  struct rv_use_class {
    // Structure describing where a rate vector can apply.
    std::string state;
    std::list<int> pos;
  };

  class raw_rate_vector {
  public:
    raw_rate_vector(std::string, rv_use_class, std::list<raw_param>);
    int ID;
    std::string name;
    rv_use_class uc;
    std::list<raw_param> rates;
    friend std::ostream& operator<<(std::ostream&, const IO::raw_rate_vector&);
  };
  
  class raw_substitution_model {
  public:
    std::string name;
    raw_substitution_model();
    void read_from_file(std::ifstream&);
    friend std::ostream& operator<<(std::ostream&, const IO::raw_substitution_model&);
    std::list<std::string> states;
    std::list<raw_rate_vector> rv_list;
    std::list<raw_param> get_parameters();
  private:
    IO::raw_param new_parameter(std::string, sol::table tbl);
    IO::raw_rate_vector new_rate_vector(std::string, sol::table, sol::table);
    void set_states(sol::table tbl);
    void add_rate_vector(raw_rate_vector rv);
  };

  raw_substitution_model* read_substitution_model(std::ifstream&);
}
#endif

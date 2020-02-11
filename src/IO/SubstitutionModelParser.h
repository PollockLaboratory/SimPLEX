#ifndef SubstitutionModelParser_h_
#define SubstitutionModelParser_h_

#include <string>
#include <iostream>
#include <list>
#include <map>

#include "sol2/sol.hpp"

namespace IO {
  enum param_type { CATEGORY, FLOAT };

  class raw_param {
  public:
    raw_param(std::string, param_type);
    int ID;
    std::string name;
    param_type t;
    friend std::ostream& operator<<(std::ostream&, const IO::raw_param&);
    virtual void read_options_table(sol::table tbl)=0;
  };

  class raw_ContinuousFloat : public raw_param {
  public:
    double init;
    double step_size;
    raw_ContinuousFloat(std::string, param_type, sol::table);
    void read_options_table(sol::table tbl);
  };

  class raw_DiscreteFloat : public raw_param {
  public:
    std::vector<float> categories;
    raw_DiscreteFloat(std::string, param_type, sol::table);
    void read_options_table(sol::table tbl);
  };

  class raw_param_wrapper {
  public:
    raw_param_wrapper(raw_param*);
    raw_param* ptr;
  };

  // RATE VECTOR

  struct rv_use_class {
    // Structure describing where a rate vector can apply.
    std::string state;
    std::list<int> pos;
  };

  class raw_rate_vector {
  public:
    raw_rate_vector(std::string, rv_use_class, std::list<raw_param*>);
    int ID;
    std::string name;
    rv_use_class uc;
    std::list<raw_param*> rates;
    friend std::ostream& operator<<(std::ostream&, const IO::raw_rate_vector&);
  };
  
  class raw_substitution_model {
  public:
    std::string name;
    raw_substitution_model();
    void read_from_file(std::ifstream&);
    friend std::ostream& operator<<(std::ostream&, const IO::raw_substitution_model&);
    std::list<std::string> states;
    std::list<std::string> ignore_states;
    std::list<raw_rate_vector> rv_list;
    std::list<raw_param*> get_parameters();
    const std::map<int, raw_param*>& get_map_parameters();
  private:
    std::map<int, raw_param*> params;
    IO::raw_param* new_parameter(std::string, std::string, sol::table tbl);
    IO::raw_rate_vector new_rate_vector(std::string, sol::table, sol::table);
    void set_states(sol::table tbl);
    void set_ignore_states(sol::table tbl);
    void add_rate_vector(raw_rate_vector rv);
  };

  raw_substitution_model* read_substitution_model(std::ifstream&);
}

#endif

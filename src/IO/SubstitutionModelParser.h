#ifndef SubstitutionModelParser_h_
#define SubstitutionModelParser_h_

#define SOL_EXCEPTIONS_SAFE_PROPAGATION

#include <string>
#include <list>
#include <map>

#include "sol/sol.hpp"
#include "SubstitutionModelParameterWrapper.h"
#include "SequencesParser.h"

namespace IO {
  struct StateData {
    enum Tag { DYNAMIC } tag;
    union {
      const RawMSA* dynamic;
    } data;

    ~StateData() {};
  };
    
  // RATE VECTOR
  struct RVScope { // RENAME ME
    // Structure describing where a rate vector applies.
    std::string domain;
    std::string state;
    std::map<std::string, std::string> secondary_state;
  };

  class raw_rate_vector {
  public:
    raw_rate_vector(std::string, RVScope, std::list<AbstractComponent*>);
    int ID;
    std::string name;
    RVScope uc;
    std::list<AbstractComponent*> rates;
    friend std::ostream& operator<<(std::ostream&, const IO::raw_rate_vector&);
  };
  
  class raw_substitution_model {
  public:
    std::string name;
    raw_substitution_model();
    void read_from_file(std::string file_name);
    
    const std::map<std::string, std::list<std::string>> get_all_states();
    const std::list<raw_rate_vector> get_rate_vector_list();

    const StateData get_state_data(std::string);

    std::map<std::string, std::list<std::string>> all_states;
    std::map<std::string, std::string> states_seqs_output_files;
    std::map<std::string, std::string> states_subs_output_files;

    std::list<raw_rate_vector> rv_list;
  private:
    void add_rate_vector(raw_rate_vector rv);

    //std::map<std::string, const IO::RawMSA*> state_data;
    std::map<std::string, StateData> state_data;

    void add_state(std::string name, sol::table states, sol::table options);
    void read_state_file(std::string domain, std::string file_name);
    void generate_uniform_data(std::string domain);
  };

  raw_substitution_model* read_substitution_model(std::string file_name);
}

#endif

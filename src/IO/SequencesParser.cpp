#include "SequencesParser.h"
#include <iostream>
#include <algorithm>
#include <exception>

#include <set>

extern IO::Files files;

std::string cleanLine(std::string line) {
  /* 
   * Cleans up the a line to remove blank spaces a line returns.
   */
  //Remove spaces and concatenates words
  line.erase(std::remove(line.begin(), line.end(), ' '), line.end());
  //This removes the carriage return in Windows files
  line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
  return(line);
}

void addSequence(IO::RawMSA& msa, std::string name, std::string seq) {
  if(name == "") {
    std::cerr << "Error: missing name in sequence alignement." << std::endl;
    exit(EXIT_FAILURE);
  }

  if(seq == "") {
    std::cerr << "Error: missing sequence in sequence alignement." << std::endl;
    exit(EXIT_FAILURE);
  }

  if(msa.n == 0) {
    msa.cols = seq.size();
  } else {
    if(seq.size() != msa.cols) {
      std::cerr << "Error: the sequence \"" << name << "\" is not the expected length." << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  msa.n++;
  msa.seqs[name] = seq;
}

namespace IO {
  // Operators
  bool operator==(const RawMSA& lhs, const RawMSA& rhs) {
    if(lhs.n != rhs.n or lhs.cols != rhs.cols) {
      return(false);
    }

    for(auto it = lhs.seqs.begin(); it != lhs.seqs.end(); ++it) {
      if((*it).second != rhs.seqs.at((*it).first)) {
	//std::cerr << (*it).second << " " << rhs.seqs.at((*it).first) << std::endl;
	return(false);
      }
    }

    return(true); 
  }

  std::ostream& operator<<(std::ostream& os, const RawMSA& msa) {
    os << "<< MSA: n:" << msa.n << " cols:" << msa.cols << " >>" << std::endl;
    for(auto it = msa.seqs.begin(); it != msa.seqs.end(); ++it) {
      os << (*it).first << " " << (*it).second << std::endl;
    }
    return(os);
  }

  RawMSA* readRawMSA(std::string file_name) {
    RawMSA* raw_msa = new RawMSA();
    raw_msa->n = 0;

    std::string name = "";
    std::string seq = "";
    std::string line; // = files.get_next_line("sequences_in");

    while(not files.end("sequences_in")) {
      line = files.get_next_line("sequences_in");
      std::string clean_line = cleanLine(line);

      //Ignore empty lines or comments.
      if(clean_line == "" or clean_line.at(0) == '#') {
	continue;
      }

      if(clean_line.at(0) == '>') {
	if (seq != "") {
	  addSequence(*raw_msa, name, seq);
	}

	name = clean_line.substr(1);
	seq = "";
      } else {
	seq += clean_line;
      }
    }
    addSequence(*raw_msa, name, seq);
    return(raw_msa);
  }

  void printRawMSA(const RawMSA& msa) {
    std::cout << "Raw Multiple Sequence Alignment:" << std::endl;
    for(auto it = msa.seqs.begin(); it != msa.seqs.end(); ++it) {
      std::cout << (*it).first << " " << (*it).second << std::endl;
    }
  }

  std::list<std::string> getRawMSANames(const RawMSA& msa) {
    std::list<std::string> names = {};
    for(auto it = msa.seqs.rbegin(); it != msa.seqs.rend(); ++it) {
      names.push_front((*it).first);
    }
    return(names);
  }

  void convertToGaps(RawMSA& msa, std::list<std::string> remove_list) {
    for(auto it = msa.seqs.begin(); it != msa.seqs.end(); ++it) {
      std::string *seq = &it->second;
      for(auto jt = seq->begin(); jt != seq->end(); ++jt) {
	if(std::find(remove_list.begin(), remove_list.end(), std::string(1, *jt)) != remove_list.end()) {
	  // Set to gap.
	  *jt = '-';
	}
      }
    }
  }

  // New structs.

  std::set<std::string> readStates(std::list<std::string> states_list) {
    std::set<std::string> states = {};
    std::set<std::string> reserved = {">", "<", ":", "-", "[", "]", ",", ";"};
 
    for(auto it = states_list.begin(); it != states_list.end(); ++it) {
      if(reserved.find(*it) != reserved.end()) {
	throw ParseException("\"" + *it + "\" is a reserved char and cannot be used as a state");
      } else {
	states.insert(*it);
      }
    }
							     
    return(states);
  }
  
  std::string freqListAsStr(std::list<StateFreq> freqs) {
    if(freqs.size() == 1) {
      return(std::string(1, freqs.front().state));
    } else {
      return("Fix me");
    }
  }

  std::string sequenceAsStr(FreqSequence seq) {
    std::string out = "";
    for(auto it = seq.begin(); it != seq.end(); ++it) {
      out += freqListAsStr(*it);
    }
    return(out);
  }

  bool operator==(const StateFreq& lhs, const StateFreq& rhs) {
    if(lhs.freq == rhs.freq and lhs.state == rhs.state) {
      return(true);
    } else {
      return(false);
    }
  }
  
  bool operator==(const std::list<StateFreq>& lhs, const std::list<StateFreq>& rhs) {
    if(lhs.size() != rhs.size()) {
      return(false);
    }

    for(auto it = lhs.begin(); it != lhs.end(); ++it) {
      auto jt = rhs.begin();
      bool match = false;
      while(match == false) {
	if(*it == *jt) {
	  match = true;
	} else {
	  ++jt;
	  if(jt == rhs.end()) {
	    return(false);
	  }
	}
      }
    }
    return(true);
  }

  bool operator==(const FreqSequence& lhs, const FreqSequence& rhs) {
    if(lhs.size() != rhs.size()) {
      return(false);
    }

    auto lit = lhs.begin();
    auto rit = rhs.begin();

    while(lit != lhs.end() and rit != rhs.end()) {
      if(*lit != *rit) {
	return(false);
      }
      lit++;
      rit++;
    }

    return(true);
  }

  bool operator==(const RawAdvMSA& lhs, const RawAdvMSA& rhs) {
    if(lhs.n != rhs.n) {
      return(false);
    }

    if(lhs.cols != rhs.cols) {
      return(false);
    }

    if(lhs.seqs.size() != rhs.seqs.size()){
      return(false);
    }

    for(auto it = lhs.seqs.begin(); it != lhs.seqs.end(); ++it) {
      auto jt = rhs.seqs.find(it->first);
      if(jt == rhs.seqs.end()) {
	return(false);
      }

      if(it->second != jt->second) {
	return(false);
      }
    }

    return(true);
  }

  ParseException::ParseException(std::string s) : std::invalid_argument(s) {
  }

  enum TokenType { NAME, POS };

  enum ParserState { SEQNAME, SEQUENCE };

  struct Token {
    std::string val;
    TokenType t;
  };

  // New Parser.
  std::set<char> separator = { '\t', '\n', ' ' };

  struct Token next_seqname(std::string::iterator &cur, std::string::iterator end, ParserState &state) {
    std::string val = "";
    while(separator.find(*cur) == separator.end()) {
      val += *cur;
      ++cur; 
    }

    struct Token t = {val, NAME};
    return(t);
  }

  struct Token next_position(std::string::iterator &cur) {
    std::string val = "";
    val += *cur;
    cur++;

    struct Token t = {val, POS};
    return(t);
  }
  
  struct Token next_token(std::string::iterator &cur, std::string::iterator end, ParserState &state) {
    struct Token t;
    if(state == SEQNAME) {
      t = next_seqname(cur, end, state);
      state = SEQUENCE;
    } else if(state == SEQUENCE) {
      t = next_position(cur);
      state = SEQNAME;
    }

    if(separator.find(*cur) != separator.end()) {
      while(separator.find(*cur) != separator.end()) {
	++cur;
      }
    } else {
      state = SEQUENCE;
    }

    return(t);
  }

  std::string validate_name(std::string token_val) {
    if(token_val[0] != '>') {
      throw ParseException("expecting sequence name, possible missing \">\"");
    }
    return(token_val.substr(1, token_val.size()));
  }

  std::list<StateFreq> validate_pos(std::string token_val) {
    // Only deals with positions that are frequency 1.0.
    StateFreq pos = {token_val.front(), 1.0};
    std::list<StateFreq> freqs = {pos};
    return(freqs);
  }

  unsigned int countPositions(std::map<std::string, FreqSequence> seqs) {
    if(seqs.empty()) {
      throw ParseException("empty file");
    }

    unsigned int count = 0;
    for(auto it = seqs.begin(); it != seqs.end(); ++it) {
      unsigned int ncols = it->second.size();
      if(ncols == 0) {
	throw ParseException("sequence for " + it->first + " is empty");
      }

      if(count == 0) {
	count = ncols;
      } else if(ncols != count) {
	throw ParseException("sequences in fasta file are not equal length");
      }
    }
    return(count);
  }
  
  RawAdvMSA parseRawAdvMSA(std::string data) {
    std::cout << "Reading file: " << data.size() << std::endl;
    RawAdvMSA msa = {0, 0, {}};
    auto loc = data.begin();
    ParserState state = SEQNAME;

    std::string seq_name = "";
    FreqSequence seq = {};
    while(loc != data.end()) {
      struct Token tok = next_token(loc, data.end(), state);
      if(tok.t == NAME) {
	if(seq_name != "") {
	  msa.n += 1;
	  // Add previous element to MSA.
	  msa.seqs[seq_name] = seq;
	  seq = {};
	}

	// New name-sequence pair. 
	seq_name = validate_name(tok.val);
      } else if (tok.t == POS) {
	std::list<StateFreq> freqs = validate_pos(tok.val);
	seq.push_back(freqs);
      }
    }

    if(seq_name != "") {
      msa.n += 1;
      msa.seqs[seq_name] = seq;
    }

    // Check sequences length.
    unsigned int ncols = countPositions(msa.seqs);
    msa.cols = ncols;

    return(msa);
  }

  // Validate MSA.
  bool validSequence(FreqSequence seq, std::set<std::string> states) {
    for(auto pos = seq.begin(); pos != seq.end(); ++pos) {
      for(auto s = pos->begin(); s != pos->end(); ++s) {
	if(states.find(std::string(1, s->state)) == states.end()) {
	  return(false);
	}
      }
    }

    return(true);
  }

  bool validateMSA(RawAdvMSA msa, std::set<std::string> states) {
    for(auto it = msa.seqs.begin(); it != msa.seqs.end(); ++it) {
      if(not validSequence(it->second, states)) {
	throw ParseException("sequence for " + it->first + " contains unrecognized state");
      }
    }
    return(true);
  }

  RawAdvMSA readRawAdvMSA(std::string data, std::set<std::string> states) {
    RawAdvMSA msa = parseRawAdvMSA(data);

    validateMSA(msa, states);

    return(msa);
  }

  void printRawAdvMSA(RawAdvMSA msa) {
   for(auto it = msa.seqs.begin(); it != msa.seqs.end(); ++it) {
     std::cout << it->first << " - " << sequenceAsStr(it->second) << std::endl;
   }
 }
}




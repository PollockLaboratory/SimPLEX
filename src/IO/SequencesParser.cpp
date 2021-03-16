#include <iostream>
#include <algorithm>
#include <exception>
#include <stdlib.h>

#include "SequencesParser.h"
#include "Files.h"

extern IO::Files files;

namespace IO {
  std::list<std::string> readStates(std::list<std::string> states_list) {
    std::list<std::string> states = {};
    std::set<std::string> reserved = {">", "<", ":", "-", "[", "]", ",", ";", "(", ")", "*"};
 
    for(auto it = states_list.begin(); it != states_list.end(); ++it) {
      if(reserved.find(*it) != reserved.end()) {
	throw ParseException("\"" + *it + "\" is a reserved char and cannot be used as a state");
      } else {
	states.push_back(*it);
      }
    }
							     
    return(states);
  }
  
  std::string freqListAsStr(std::list<StateFreq> freqs) {
    if(freqs.size() == 1) {
      return(std::string(1, freqs.front().state));
    } else {
      float max_freq = 0.0;
      char state = '*';
      for(auto it = freqs.begin(); it != freqs.end(); ++it) {
	if(it->freq > max_freq) {
	  max_freq = it->freq;
	  state = it->state;
	}
      }
      return(std::string(1, state));
    }
  }

  bool FreqList_gap(std::list<StateFreq> freqs) {
    if(freqs.empty()) {
      throw ParseException("empty list of state frequencies.");
    } else if(freqs.size() >= 2) {
      for(auto it = freqs.begin(); it != freqs.end(); ++it) {
	if(it->state == '-') {
	  throw ParseException("position specifies gap at frequency less that 1.0.");
	}
      }
    } else {
      if(freqs.front().state == '-') {
	return(true);
      }
    }
    return(false);
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

  bool operator==(const RawMSA& lhs, const RawMSA& rhs) {
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

  enum TokenType { NAME, POS, POSEXT };

  enum ParserState { SEQNAME, SEQUENCE, SEQUENCEEXT };

  struct Token {
    std::string val;
    TokenType t;
    float freq;
  };

  // New Parser.
  std::set<char> separator = { '\t', '\n', ' ' };

  struct Token next_seqname(std::string::iterator &cur, std::string::iterator end, ParserState &state) {
    std::string val = "";
    while(separator.find(*cur) == separator.end()) {
      val += *cur;
      ++cur; 
    }

    struct Token t = {val, NAME, 0.0};
    return(t);
  }

  struct Token next_position(std::string::iterator &cur) {
    std::string val = "";

    if(*cur == '[') {
      cur++;
      while(*cur != ']' and *cur != ',') {
	val += *cur;
	cur++;
      }

      size_t split = val.find(':');
      std::string state = val.substr(0,split);
      float freq = std::stod(val.substr(split+1));

      struct Token t;
      if(*cur == ',') {
	t = {state, POSEXT, freq};
      } else if (*cur == ']') {
	t = {state, POS, freq};
      }

      cur++;

      return(t);
    } else {
      val += *cur;
      cur++;
      
      struct Token t = {val, POS, 1.0};
      return(t);
    }
  }

  struct Token next_position_ext(std::string::iterator &cur) {
    std::string val = "";
    while(*cur != ']' and *cur != ',') {
      val += *cur;
      cur++;
    }

    size_t split = val.find(':');
    std::string state = val.substr(0,split);
    float freq = std::stod(val.substr(split+1));

    struct Token t;
    if(*cur == ',') {
      t = {state, POSEXT, freq};
    } else if (*cur == ']') {
      t = {state, POS, freq};
    }

    cur++;
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
    } else if(state == SEQUENCEEXT) {
      t = next_position_ext(cur);
      if(t.t == POS) {
	state = SEQNAME;
      }
    }

    // Do complex logic here.
    if(t.t == POSEXT) {
      state = SEQUENCEEXT;
    } else if(separator.find(*cur) != separator.end()) {
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
  
  RawMSA parseRawAdvMSA(std::string data) {
    RawMSA msa = {0, 0, {}};
    auto loc = data.begin();
    ParserState state = SEQNAME;

    std::string seq_name = "";
    FreqSequence seq = {};
    std::list<StateFreq> freqs = {};
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
      } else if(tok.t == POS) {
	// Push onto seq.
	StateFreq pos = {tok.val.front(), tok.freq};
	freqs.push_back(pos);
	seq.push_back(freqs);
	freqs = {};
      } else if(tok.t == POSEXT) {
	StateFreq pos = {tok.val.front(), tok.freq};
	freqs.push_back(pos);
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
	  // Need something better here.
	  if(s->state != '-') {
	    return(false);
	  }
	}
      }
    }

    return(true);
  }

  bool validateMSA(RawMSA msa, std::set<std::string> states) {
    for(auto it = msa.seqs.begin(); it != msa.seqs.end(); ++it) {
      if(not validSequence(it->second, states)) {
	throw ParseException("sequence for " + it->first + " contains unrecognized state");
      }
    }
    return(true);
  }

  bool checkEmptyColumns(RawMSA msa) {
    std::vector<bool> empty_cols(msa.cols, true);
    for(auto seq = msa.seqs.begin(); seq != msa.seqs.end(); ++seq) {
      unsigned int i = 0;
      for(auto pos = seq->second.begin(); pos != seq->second.end(); ++pos) {
	if(not FreqList_gap(*pos)) {
	  empty_cols[i] = false;
	}
	i++;
      }
    }
    for(auto it = empty_cols.begin(); it != empty_cols.end(); ++it) {
      if(*it == true) {
	throw ParseException("empty column (columns containing only gaps)");
      }
    }
  }

  RawMSA readRawAdvMSA(std::string data, std::list<std::string> states) {
    RawMSA msa = parseRawAdvMSA(data);

    std::set<std::string> states_set(states.begin(), states.end());
    validateMSA(msa, states_set);
    checkEmptyColumns(msa);

    return(msa);
  }

  void printRawAdvMSA(RawMSA msa) {
   for(auto it = msa.seqs.begin(); it != msa.seqs.end(); ++it) {
     std::cout << it->first << " - " << sequenceAsStr(it->second) << std::endl;
   }
 }
}




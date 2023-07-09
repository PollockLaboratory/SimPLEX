#include <iostream>
#include <algorithm>
#include <iterator>
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
  
  std::string freqListAsStr_highestFreq(std::list<StateFreq> freqs) {
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
 
  std::string freqListAsStr(std::list<StateFreq> freqs) {
    if(freqs.size() == 1) {
      return(std::string(1, freqs.front().state));
    } else {
      std::string freq_str = "[";
      for(auto it = freqs.begin(); it != freqs.end(); ++it) {
        freq_str += std::string(1, it->state) + ":" + std::to_string(it->freq);
        if(std::next(it) != freqs.end()) {
          freq_str += ",";
        }
      }
      freq_str += "]";
      
      return(freq_str);
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

  std::string sequenceAsStr_highestFreq(FreqSequence seq) {
    // Prints a string representation of the FreqSequences, but picks the highest frequences
    // state to represent the position.
    std::string out = "";
    for(auto it = seq.begin(); it != seq.end(); ++it) {
      out += freqListAsStr_highestFreq(*it);
    }
    return(out);
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
    TokenType tag;
    float freq;
  };

  // New Parser.
  std::set<char> separator = { '\t', '\n', ' ' };

  struct Token next_seqname(std::string::iterator &cur) {
    std::string val = "";
    while(separator.find(*cur) == separator.end()) {
      val += *cur;
      ++cur; 
    }

    return(Token({val, NAME, 0.0}));
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
      
      return(Token({val, POS, 1.0}));
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
  
  struct Token next_token(std::string::iterator &cur, ParserState &state) {
    struct Token next_token;
    if(state == SEQNAME) {
      next_token = next_seqname(cur);
      state = SEQUENCE;
    } else if(state == SEQUENCE) {
      next_token = next_position(cur);
      state = SEQNAME;
    } else if(state == SEQUENCEEXT) {
      next_token = next_position_ext(cur);
      if(next_token.tag == POS) {
        state = SEQNAME;
      }
    }

    // Do complex logic here.
    if(next_token.tag == POSEXT) {
      state = SEQUENCEEXT;
    } else if(separator.find(*cur) != separator.end()) {
      while(separator.find(*cur) != separator.end()) {
        ++cur;
      }
    } else {
      state = SEQUENCE;
    }

    return(next_token);
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
  
  RawMSA parseRawMSA(std::string data) {
    // Parses MSA with frquency infomation.
    // A = state 'A' with frequency of 1.0.
    // [A:1.0] = state 'A' with frequency of 1.0.
    // [A:0.5,B:0.5] = state 'A' with frequncy of 0.5 and state 'B' with frequency of 0.5.
    RawMSA msa = {0, 0, {}};
    auto loc = data.begin();
    ParserState state = SEQNAME;

    std::string seq_name = "";
    FreqSequence seq = {};
    std::list<StateFreq> freqs = {};
    while(loc != data.end()) {
      struct Token tok = next_token(loc, state);
      if(tok.tag == NAME) {
        if(seq_name != "") {
          msa.n += 1;
          // Add previous element to MSA.
          msa.seqs[seq_name] = seq;
          seq = {};
        }

        // New name-sequence pair. 
        seq_name = validate_name(tok.val);
      } else if(tok.tag == POS) {
        // Push onto seq.
        StateFreq pos = {tok.val.front(), tok.freq};
        freqs.push_back(pos);
        seq.push_back(freqs);
        freqs = {};
      } else if(tok.tag == POSEXT) {
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

  void validateMSA(RawMSA msa, std::set<std::string> states) {
    for(auto it = msa.seqs.begin(); it != msa.seqs.end(); ++it) {
      if(not validSequence(it->second, states)) {
        throw ParseException("sequence for " + it->first + " contains unrecognized state");
      }
    }
  }

  void checkEmptyColumns(RawMSA msa) {
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

    // throw error if column is only gaps.
    for(auto it = empty_cols.begin(); it != empty_cols.end(); ++it) {
      if(*it) {
        throw ParseException("empty column (columns containing only gaps)");
      }
    }
  }

  RawMSA readRawMSA(std::string data, std::list<std::string> states) {
    // This parses a .fasta or efasta file (extended fasta).
    // .efasta are a superset of fasta file format.
    //
    // In fasta all positions are represented as single characters, representing
    // a state with a frequency of 1.0, efasta can include uncertains states.
    // To represent an uncertain state use [state:freq, ... ].
    // For example; [A:0.5, C:0.5] where a position has 0.5 frequency of either A or C state.
    // Note: A and [A:1.0] are equivilent.
    // Note: whitespace is ignores within square brackets.
    //
    // example of name and sequence pair.
    //>taxaA
    // ATC[A:0.5,C:0.5]
    //
    RawMSA msa = parseRawMSA(data);

    std::set<std::string> states_set(states.begin(), states.end());
    validateMSA(msa, states_set);
    checkEmptyColumns(msa);

    return(msa);
  }

  RawMSA createUniformPrior(std::list<std::string> states, RawMSA template_msa) {
    RawMSA msa = {};

    // Construct uniform frequences for single position.
    float uniform_freq = 1.0 / states.size();
    std::list<StateFreq> uniform_pos = {};
    for(auto state_it = states.begin(); state_it != states.end(); ++state_it) {
      StateFreq sf = {state_it->front(), uniform_freq};
      uniform_pos.push_back(sf);
    }

    std::list<StateFreq> gap_pos = {{'-', 1.0}};

    for(auto it = template_msa.seqs.begin(); it != template_msa.seqs.end(); ++it) {
      std::string seq_name = it->first;
      FreqSequence seq = it->second;

      msa.seqs[seq_name] = {};
      // Loop through positions of template MSA, if gap insert gap else add uniform prior.
      for(auto pos_it = seq.begin(); pos_it != seq.end(); ++pos_it) {
	if(FreqList_gap(*pos_it)) {
	  msa.seqs[seq_name].push_back(gap_pos);
	} else {
	  msa.seqs[seq_name].push_back(uniform_pos);
	}
      }
    }

    return(msa);
  }

  void printRawAdvMSA(RawMSA msa) {
   for(auto it = msa.seqs.begin(); it != msa.seqs.end(); ++it) {
     std::cout << it->first << " - ";
     std::string seq_str = sequenceAsStr(it->second);
     if(seq_str.empty()) {
       std::cout << "EMPTY_SEQUENCE" << std::endl;
     } else {
       std::cout << seq_str << std::endl;
     }
   }
 }
}

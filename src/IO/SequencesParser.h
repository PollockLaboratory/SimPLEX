#ifndef SequencesParser_h_
#define SequencesParser_h_

#include <map>
#include <list>
#include <set>

namespace IO {
  class ParseException : public std::invalid_argument {
  public:
    ParseException(std::string);
  };

  std::list<std::string> readStates(std::list<std::string> states);
  
  typedef struct {
    char state;
    float freq;
  } StateFreq;

  std::string freqListAsStr_highestFreq(std::list<StateFreq> freqs);
  std::string freqListAsStr(std::list<StateFreq> freqs);
  bool FreqList_gap(std::list<StateFreq> freqs);

  // Type for representing sequences where each position is contains state frequencies.
  typedef std::list<std::list<StateFreq>> FreqSequence;
  bool operator==(const FreqSequence& lhs, const FreqSequence& rhs);
  std::string sequenceAsStr_highestFreq(FreqSequence seq);
  std::string sequenceAsStr(FreqSequence seq);
  
  struct RawMSA {
    unsigned int n = 0;
    unsigned int cols = 0;
    std::map<std::string, FreqSequence> seqs;
  };

  RawMSA parseRawMSA(std::string data);
  RawMSA readRawMSA(std::string data, std::list<std::string> states);
  RawMSA createUniformPrior(std::list<std::string> states, RawMSA template_msa);

  // Utils
  bool operator==(const RawMSA& lhs, const RawMSA& rhs);
  void printRawAdvMSA(RawMSA msa);
}

#endif

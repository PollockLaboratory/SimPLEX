#ifndef SequencesParser_h_
#define SequencesParser_h_

#include <map>
#include <list>
#include <set>

namespace IO {
  struct RawMSA {
    std::map<std::string, std::string> seqs;
    unsigned int n = 0;
    unsigned int cols = 0;
  };

  // Operators
  bool operator==(const RawMSA& lhs, const RawMSA& rhs);
  std::ostream& operator<<(std::ostream& os, const RawMSA& msa);

  // Utils
  RawMSA* readRawMSA(std::string file_name);
  void printRawMSA(const RawMSA& msa);
  std::list<std::string> getRawMSANames(const RawMSA& msa);
  void convertToGaps(RawMSA& msa, std::list<std::string> remove_list);

  // New functions and struct.

  class ParseException : public std::invalid_argument {
  public:
    ParseException(std::string);
  };

  std::list<std::string> readStates(std::list<std::string> states);
  
  typedef struct {
    char state;
    float freq;
  } StateFreq;

  std::string freqListAsStr(std::list<StateFreq> freqs);
  bool FreqList_gap(std::list<StateFreq> freqs);

  // Type for representing sequences where each position is contains state frequencies.
  typedef std::list<std::list<StateFreq>> FreqSequence;
  bool operator==(const FreqSequence& lhs, const FreqSequence& rhs);
  std::string sequenceAsStr(FreqSequence seq);
  
  struct RawAdvMSA {
    unsigned int n = 0;
    unsigned int cols = 0;
    std::map<std::string, FreqSequence> seqs;
  };

  RawAdvMSA parseRawAdvMSA(std::string data);
  RawAdvMSA readRawAdvMSA(std::string data, std::list<std::string> states);

  // Utils
  bool operator==(const RawAdvMSA& lhs, const RawAdvMSA& rhs);
  void printRawAdvMSA(RawAdvMSA msa);
}

#endif

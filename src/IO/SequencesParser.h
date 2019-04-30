#ifndef SequencesParser_h_
#define SequencesParser_h_

#include <map>
#include <list>
#include <algorithm>
#include "Files.h"

namespace IO {
  struct RawMSA {
    std::map<std::string, std::string> seqs;
    unsigned int n = 0;
    unsigned int cols = 0;
  };

  RawMSA* readRawMSA(std::ifstream &sequences_stream);
  void printRawMSA(const RawMSA& msa);
  std::list<std::string> getRawMSANames(const RawMSA& msa);
}

#endif
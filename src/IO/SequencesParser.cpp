#include "SequencesParser.h"
#include <iostream>
#include <algorithm>

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
    for(auto it = msa.seqs.begin(); it != msa.seqs.end(); ++it) {
      names.push_front((*it).first);
    }
    return(names);
  }

  void convertToGaps(RawMSA& msa, std::list<std::string> remove_list) {
    for(auto it = msa.seqs.begin(); it != msa.seqs.end(); ++it) {
      std::string *seq = &it->second;
      for(auto jt = seq->begin(); jt != seq->end(); ++jt) {
	if(std::find(remove_list.begin(), remove_list.end(), std::string(1, *jt)) != remove_list.end()) {
	  *jt = '-';
	}
      }
    }
  }
}

// Utils


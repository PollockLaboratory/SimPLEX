#ifndef ObservedData_h_
#define ObservedData_h_

#include <iostream>
#include <fstream>

#include <string>
#include <vector>
#include <set> // for columns_with_gaps
#include <map> // for taxa_names_to_sequences
#include <list>
#include <algorithm>

#include "Sequence.h"
#include "Trees/TreeParser.h"
#include "SubstitutionModels/SubstitutionModel.h"

using namespace std;

class Data {
 public:
  SequenceAlignment* MSA;
  list<SequenceAlignment*> MSA_list;
  IO::RawTreeNode* raw_tree;

  Data();
  ~Data();

  void Initialize(const States*);
 private:
  set<int> columns_with_gaps;
  vector<int> columns_without_gaps;

  string cleanLine(string);
  list<string> readFastaFile(ifstream &sequences_file);
  list<list<string>> readCompoundFastaFile(ifstream &sequences_file);

  void validateInputData(list<SequenceAlignment*> MSA_list, IO::RawTreeNode* raw_tree);
  bool matchNodeNames(list<string> names1, list<string> names2);

  SequenceAlignment* ReadSequences(list<string> fasta_lines, const States*);

  IO::RawTreeNode* ReadTree();
};

#endif


#ifndef Data_h_
#define Data_h_

#include <iostream>
#include <fstream>

#include <string>
#include <vector>
#include <set> // for columns_with_gaps
#include <map> // for taxa_names_to_sequences
#include <list>
#include <algorithm>

#include "ModelParts/Sequence.h"
#include "IO/TreeParser.h"
#include "ModelParts/SubstitutionModels/SubstitutionModel.h"
#include "IO/SubstitutionModelParser.h"
#include "IO/SequencesParser.h"

using namespace std;

class Data {
 public:
  SubstitutionModel* sm;
  SequenceAlignment* MSA;
  list<SequenceAlignment*> MSA_list;
  IO::RawMSA* raw_msa;
  IO::RawTreeNode* raw_tree;
  IO::raw_substitution_model* raw_sm;

  Data();
  ~Data();

  void Initialize();
 private:
  set<int> columns_with_gaps;
  vector<int> columns_without_gaps;

  void validateInputData(const IO::RawMSA* raw_msa, const IO::RawTreeNode* raw_tree);
  void matchNodeNames(list<string> names1, list<string> names2);
  
  IO::RawTreeNode* ReadTree();
  IO::RawMSA* ReadMSA();
  IO::raw_substitution_model* ReadSubstitutionModel(const IO::RawMSA*, const IO::RawTreeNode*);
};

#endif


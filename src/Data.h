#ifndef Data_h_
#define Data_h_

#include <string>
#include <vector>
#include <set> // for columns_with_gaps
#include <list>
#include <algorithm>

#include "IO/TreeParser.h"
#include "IO/SubstitutionModelParser.h"
#include "IO/SequencesParser.h"

class Data {
public:
  IO::RawTreeNode* raw_tree;
  IO::raw_substitution_model* raw_sm;

  Data();
  ~Data();

  void Initialize();
  void Uninitialize();
private:
  IO::RawTreeNode* ReadTree();
  IO::raw_substitution_model* ReadSubstitutionModel();
};

#endif


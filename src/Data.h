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
  IO::RawMSA* raw_msa;
  IO::RawTreeNode* raw_tree;
  IO::raw_substitution_model* raw_sm;

  Data();
  ~Data();

  void Initialize();
  void Uninitialize();
private:
  std::set<int> columns_with_gaps;
  std::vector<int> columns_without_gaps;

  void validateInputData(const IO::RawMSA* raw_msa, const IO::RawTreeNode* raw_tree);
  void matchNodeNames(std::list<std::string> names1, std::list<std::string> names2);
  
  IO::RawTreeNode* ReadTree();
  IO::RawMSA* ReadMSA();
  IO::raw_substitution_model* ReadSubstitutionModel();
};

#endif


#include "Data.h"

#include <algorithm>
#include <map>
#include <list>

#include "Environment.h"
#include "IO/Files.h"

extern Environment env;
extern IO::Files files;

static const int gap_indicator = -1;

Data::Data() {
//	cout << "Data constructor" << endl;
}

Data::~Data() {
//	cout << "Data destructor" << endl;
}

void Data::Initialize() {
  std::cout << "Loading data:" << std::endl;
  raw_msa = ReadMSA();
  raw_tree = ReadTree();
  raw_sm = ReadSubstitutionModel(raw_msa, raw_tree);

  std::cout << "States: " << std::endl;
  for(auto it = raw_sm->states.begin(); it != raw_sm->states.end(); ++it) {
    std::cout << *it << " ";
  }
  std::cout << std::endl;

  std::cout << "Ignore states: " << std::endl;
  for(auto it = raw_sm->ignore_states.begin(); it != raw_sm->ignore_states.end(); ++it) {
    std::cout << *it << " ";
  }
  std::cout << std::endl;

  IO::convertToGaps(*raw_msa, raw_sm->ignore_states);
  
  validateInputData(raw_msa, raw_tree);
  std::cout << "Data successfully read." << std::endl;
}

IO::RawTreeNode* Data::ReadTree() {
  std::string treefile = env.get<std::string>("DATA.tree_file");
  files.add_file("tree_input", treefile, IOtype::INPUT);

  std::cout << "Reading tree from:\t" << files.get_file_info("tree_input") << std::endl;

  ifstream tree_in = files.get_ifstream("tree_input");

  string tree_string;
  getline(tree_in, tree_string);

  IO::RawTreeNode* raw_tree = IO::parseTree(tree_string);

  return(raw_tree);
}

IO::RawMSA* Data::ReadMSA() {
  files.add_file("sequences_in", env.get<std::string>("DATA.sequences_file"), IOtype::INPUT);
  ifstream sequences_in = files.get_ifstream("sequences_in");

  std::cout << "Reading sequences from:\t" << files.get_file_info("sequences_in") << std::endl;

  IO::RawMSA* raw_msa = IO::readRawMSA(sequences_in);

  return(raw_msa);  
}

IO::raw_substitution_model* Data::ReadSubstitutionModel(const IO::RawMSA* raw_msa, const IO::RawTreeNode* raw_tree) {
  files.add_file("lua_model", env.get<std::string>("DATA.substitution_model_file"), IOtype::INPUT);
  std::ifstream lua_sm_in = files.get_ifstream("lua_model");

  std::cout << "Reading Substitution model from:\t" << files.get_file_info("lua_model") << std::endl;

  IO::raw_substitution_model* raw_sm = IO::read_substitution_model(lua_sm_in);
  return(raw_sm);
}

// Validate data.

void Data::matchNodeNames(list<string> tree_nodes, list<string> msa_names) {
  for(auto name = tree_nodes.begin(); name != tree_nodes.end(); ++name) {
    auto it = std::find(msa_names.begin(), msa_names.end(), *name);
    if(it == msa_names.end()) {
      std::cerr << "Error: missing sequence for " << *name << " in MSA." << std::endl;
      exit(EXIT_FAILURE);
    }
    msa_names.erase(it);
  }
  if(not msa_names.empty()) {
    std::cerr << "Error: " << msa_names.size() << " redundant sequences in the alignment." << std::endl;
    exit(EXIT_FAILURE);
  }
}

void Data::validateInputData(const IO::RawMSA* raw_msa, const IO::RawTreeNode* raw_tree) {
  /* 
   * Check that the node names match between the tree and MSAs.
   */
  std::list<std::string> tree_names = IO::getRawTreeNodeTipNames(raw_tree);
  std::list<std::string> MSA_names = IO::getRawMSANames(*raw_msa);

  matchNodeNames(tree_names, MSA_names);
}

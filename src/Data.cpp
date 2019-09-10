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
  raw_msa = ReadMSA();
  raw_tree = ReadTree();
  raw_sm = ReadSubstitutionModel(raw_msa, raw_tree);

  validateInputData(raw_msa, raw_tree);  
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

bool Data::matchNodeNames(list<string> names1, list<string> names2) {
  auto n1 = names1.begin();
  auto n2 = names2.begin();
  std::cout << names1.size() << " " << names2.size() << std::endl;
  for(int i = 0; i < names1.size(); i++) {
    std::cout << *n1 << " " << *n2 << std::endl;
    ++n1;
    ++n2;
    if(n1 == names1.end() or n2 == names2.end()) {
      break;
    }
  }
  
  for(auto name_it = names1.begin(); name_it != names1.end(); ++name_it) {
    auto it = std::find(names2.begin(), names2.end(), *name_it);
    if(it == names2.end()) {
      return(false);
    }
    names2.erase(it);
  }
  return(names2.empty());
}

void Data::validateInputData(const IO::RawMSA* raw_msa, const IO::RawTreeNode* raw_tree) {
  /* 
   * Check that the node names match between the tree and MSAs.
   */
  // Need to impliment for non-ancestral sequences case.
  std::list<std::string> tree_names = IO::getRawTreeNodeNames(raw_tree);
  std::list<std::string> MSA_names = IO::getRawMSANames(*raw_msa);
   //for(auto it = MSA_list.begin(); it != MSA_list.end(); ++it) {
   //MSA_names = (*it)->getNodeNames();
  //   if(matchNodeNames(tree_names, MSA_names) == false) {
  //   std::cerr << "Error: Node names in MSA and tree do not match." << std::endl;
  //   exit(EXIT_FAILURE);
  //}
}

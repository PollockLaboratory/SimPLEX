#include "Data.h"

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
  raw_sm = ReadSubstitutionModel();

  raw_msa = ReadMSA();
  raw_tree = ReadTree();

  // Print States.
  std::set<std::string> states = raw_sm->get_states();
  std::cout << "States: [ ";
  for(auto it = states.begin(); it != states.end(); ++it) {
    std::cout << *it << " ";
  }
  std::cout << "]" << std::endl;

  std::list<std::string> ignore_states = raw_sm->get_ignore_states();
  std::cout << "Ignore states: [ ";
  for(auto it = ignore_states.begin(); it != ignore_states.end(); ++it) {
    std::cout << *it << " ";
  }
  std::cout << "]" << std::endl;

  // Validation.
  IO::convertToGaps(*raw_msa, ignore_states);
  
  validateInputData(raw_msa, raw_tree);
  std::cout << "Data successfully read." << std::endl;
}

void Data::Uninitialize() {
  // Raw imput data should no longer be used during MCMC. Therefore delete data to save space.
  delete raw_msa;
  IO::deleteTree(raw_tree);
}

IO::RawTreeNode* Data::ReadTree() {
  std::string treefile = env.get<std::string>("INPUT.tree_file");
  files.add_file("tree_input", treefile, IOtype::INPUT);

  std::cout << "Reading tree from:\t" << files.get_info("tree_input").file_name << std::endl;

  std::string tree_string = files.get_next_line("tree_input");
  IO::RawTreeNode* raw_tree = IO::parseTree(tree_string);

  return(raw_tree);
}

IO::RawMSA* Data::ReadMSA() {
  files.add_file("sequences_in", env.get<std::string>("INPUT.sequences_file"), IOtype::INPUT);

  std::cout << "Reading sequences from:\t" << files.get_info("sequences_in").file_name << std::endl;

  IO::RawMSA* raw_msa = IO::readRawMSA("sequences_in");

  return(raw_msa);  
}

IO::raw_substitution_model* Data::ReadSubstitutionModel() {
  files.add_file("lua_model", env.get<std::string>("MODEL.script_file"), IOtype::INPUT);
  std::cout << "Reading Substitution model from:\t" << files.get_info("lua_model").file_name << std::endl;

  IO::raw_substitution_model* raw_sm = IO::read_substitution_model("lua_model");
  return(raw_sm);
}

// Validate data.

void Data::matchNodeNames(std::list<std::string> tree_nodes, std::list<string> msa_names) {
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

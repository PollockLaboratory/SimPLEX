#include "Data.h"

#include <map>

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

  //raw_msa = ReadMSA();
  raw_tree = ReadTree();

  std::cout << "Data successfully read." << std::endl;
}

void Data::Uninitialize() {
  // Raw imput data should no longer be used during MCMC. Therefore delete data to save space.
  //delete raw_msa;
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

IO::raw_substitution_model* Data::ReadSubstitutionModel() {
  files.add_file("lua_model", env.get<std::string>("MODEL.script_file"), IOtype::INPUT);
  std::cout << "Reading Substitution model from:\t" << files.get_info("lua_model").file_name << std::endl;

  IO::raw_substitution_model* raw_sm = IO::read_substitution_model("lua_model");
  return(raw_sm);
}

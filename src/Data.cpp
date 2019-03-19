#include "Data.h"

#include <algorithm>
#include <map>
#include <list>

#include "Environment.h"
#include "IO.h"

extern Environment env;
extern IO::Files files;

static const int gap_indicator = -1;

Data::Data() {
//	cout << "Data constructor" << endl;
}

Data::~Data() {
//	cout << "Data destructor" << endl;
}

void Data::Initialize(const States* states) {
  files.add_file("sequences_in", env.get<std::string>("DATA.sequences_file"), IOtype::INPUT);
  ifstream sequences_in = files.get_ifstream("sequences_in");

  if(env.ancestral_sequences == false) {
    // Fasta file only contains terminal sequences.
    list<string> fasta_lines = readFastaFile(sequences_in);
    MSA = ReadSequences(fasta_lines, states);
    MSA->Initialize(&MSA_list);
  } else {
    //Compound fasta with internal sequences - when ancestral states are already calculated.
    list<list<string>> fasta_blocks = readCompoundFastaFile(sequences_in);
 
    for(auto it = fasta_blocks.begin(); it != fasta_blocks.end(); ++it) {
      SequenceAlignment* msa = ReadSequences(*it, states);
      msa->Initialize(&MSA_list);
      MSA_list.push_back(msa);
    }

    // The MSA_list should be checked here. That all the MSAs have the same node names.
    // MSA = new SequenceAlignment(*(MSA_list.front()));
    std::cout << "MSA_list length: " << MSA_list.size() << std::endl;
    MSA = new SequenceAlignment(*MSA_list.front());
    // MSA->Initialize(&MSA_list);
  }

  raw_tree = ReadTree();

  validateInputData(MSA_list, raw_tree);
}

// Sequences
string Data::cleanLine(string line) {
  /* 
   * Cleans up the a line to remove blank spaces a line returns.
   */
  //Remove spaces and concatenate words
  line.erase(std::remove(line.begin(), line.end(), ' '), line.end());
  //This removes the carriage return in Windows files
  line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
  return(line);
}

list<string> Data::readFastaFile(ifstream &sequences_file) {
  list<string> f;
  string line;
  while(sequences_file.good()) {
    getline(sequences_file, line);
    f.push_back(cleanLine(line));
  }
  return(f);
}

list<list<string>> Data::readCompoundFastaFile(ifstream &sequences_file) {
  list<list<string>> lf;
  list<string> f;
  string line;

  while(sequences_file.good()) {
    getline(sequences_file, line);
    line = cleanLine(line);

    // Just ignore empty lines.
    if(line != "") {
      if(line.at(0) == '#') {
	// Next fasta block.
	if(!f.empty()) {
	  lf.push_back(f);
	  f.clear();
	}
      }
      f.push_back(line);
    }
  }
  // Add last block.
  lf.push_back(f);
  return(lf);
}

SequenceAlignment* Data::ReadSequences(list<string> fasta_lines, const States* states) {
  /*
   * Given a list of lines that represent a fasta file, returns a pointer to a MSA object.
   */
  
  MSA = new SequenceAlignment(states);

  string line;
  string sequence = "";
  string name;

  for(auto it = fasta_lines.begin(); it != fasta_lines.end(); ++it) {
    line = *it;
    if(line == "" or line.at(0) == '#') {
      continue;
    }
    if(line.at(0) == '>') {
      if (sequence != "") {
	MSA->add(name, sequence);
      }
      name = line.substr(1);
      sequence = "";
    } else {
      sequence += line;
    }
  }
  // Add final sequence
  MSA->add(name, sequence);
  return(MSA);
}

// Trees

IO::RawTreeNode* Data::ReadTree() {
  std::string treefile = env.get<std::string>("DATA.tree_file");
  files.add_file("tree_input", treefile, IOtype::INPUT);

  std::cout << "Reading tree from:\t" << files.get_file_path("tree_input") << std::endl;

  ifstream tree_in = files.get_ifstream("tree_input");

  string tree_string;
  getline(tree_in, tree_string);

  IO::RawTreeNode* raw_tree = IO::parseTree(tree_string);

  return(raw_tree);
}

// Validate data
bool Data::matchNodeNames(list<string> names1, list<string> names2) {
  for(auto name_it = names1.begin(); name_it != names1.end(); ++name_it) {
    auto it = std::find(names2.begin(), names2.end(), *name_it);
    if(it == names2.end()) {
      return(false);
    }
    names2.erase(it);
  }
  return(names2.empty());
}

void Data::validateInputData(list<SequenceAlignment*> MSA_list, IO::RawTreeNode* raw_tree) {
  /* 
   * Check that the node names match between the tree and MSAs.
   */
  // Need to impliment for non-ancestral sequences case.
  if(env.ancestral_sequences) {
    list<string> tree_names = IO::getRawTreeNodeNames(raw_tree);
    list<string> MSA_names;
    for(auto it = MSA_list.begin(); it != MSA_list.end(); ++it) {
      MSA_names = (*it)->getNodeNames();
      if(matchNodeNames(tree_names, MSA_names) == false) {
	std::cerr << "Error: Node names in Multiple Sequence Alignments and tree do not match." << std::endl;
	exit(EXIT_FAILURE);
      }
    }
  }
}

#include "Data.h"

#include <algorithm>
#include <map>

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

void Data::Initialize() {
	MSA = ReadSequences();
	MSA->print();
	MSA->Initialize();

	raw_tree = ReadTree();
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
	return line;
}

SequenceAlignment* Data::ReadSequences() {
	MSA = new SequenceAlignment();

	files.add_file("sequences_in", env.get("sequences_file"), IOtype::INPUT);
	ifstream sequences_in = files.get_ifstream("sequences_in");

	vector<int> encoded_sequence;
	string line;
	string sequence = "";
	string name;

	while (sequences_in.good()) {
		getline(sequences_in, line);
		line = cleanLine(line);

		if (line == "")
			continue;

		if (line.at(0) == '>') {
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
	std::string treefile = env.get("tree_file"); 	
	files.add_file("tree_input", treefile, IOtype::INPUT);
	ifstream tree_in = files.get_ifstream("tree_input");

	string tree_string;
	getline(tree_in, tree_string);

	IO::RawTreeNode* raw_tree = IO::parseTree(tree_string);

	return(raw_tree);

}

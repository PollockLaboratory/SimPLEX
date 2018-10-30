#include "TreeParser.h"
#include "Environment.h"
#include "IO.h"

#include <algorithm>

extern Environment env;
extern IO::Files files;

inline std::string& IO::cleanTreeString(std::string &tree_string) {
	/*
	 * Removes whitespace and trailing ';'.
	 */
	tree_string.erase(remove_if(tree_string.begin(), tree_string.end(), isspace), tree_string.end());
	if(tree_string.back() == ';') {
		tree_string.erase(tree_string.end()-1, tree_string.end());
	}
	return(tree_string);
}

std::pair<std::string, std::string> IO::splitBranchString(std::string branch_string) {
	std::pair<std::string, std::string> b;
	int tree_depth = 0;
	for (int position = 0; position < branch_string.size(); position++) {
		char character = branch_string.at(position);
		if (character == '(')
			tree_depth++;
		else if (character == ')')
			tree_depth--;
		else if (character == ',' and tree_depth == 0) {
			b.first = branch_string.substr(0, position);
			b.second = branch_string.substr(position + 1, branch_string.size() - position + 1);
		}
	}
	return(b);
}

node_data IO::deconstructNodeString(std::string node_string) {
	static int ID = 0;
	std::cout << node_string << std::endl;
	int last_colon_position = node_string.find_last_of(':');
	int last_parens_position = node_string.find_last_of(')');

	float distance = atof(node_string.substr(last_colon_position + 1).c_str());	

	std::string name;
	std::pair<std::string, std::string> branch_strings;
	if (last_parens_position == string::npos) {
		// Tip node.
		name = node_string.substr(0, last_colon_position);
		branch_strings = std::make_pair("", "");
	} else {
		int len = last_colon_position - last_parens_position - 1;
		name = node_string.substr(last_parens_position + 1, len);
		if(name == "") {
			name = "BNode" + std::to_string(ID);
			ID++;
		}
		std::cout << "Name: " << name << std::endl;
		branch_strings = splitBranchString(node_string.substr(1, last_parens_position - 1));
	}
	std::cout << "Distance: " << distance << " Name: " << name << std::endl;

	node_data n = {name, distance, branch_strings.first, branch_strings.second};
	return(n);
}

IO::RawTreeNode* IO::parseRawTreeNode(std::string node_string, RawTreeNode* up) {
	node_data n = deconstructNodeString(node_string);

	IO::RawTreeNode* t = new RawTreeNode;
	IO::RawTreeNode* left;
	IO::RawTreeNode* right;

	left = (n.left != "") ? parseRawTreeNode(n.left, t) : 0;
	right = (n.right != "") ? parseRawTreeNode(n.right, t) : 0;

	*t = {n.name, n.distance, up, left, right};
	return(t);
}

IO::RawTreeNode* IO::parseTree(std::string tree_string) {
	std::cout << "Reading tree: " << std::endl;
	tree_string = cleanTreeString(tree_string);
	IO::RawTreeNode* root = new RawTreeNode;
	*root = {"Root", 0, NULL, NULL, NULL};
	IO::RawTreeNode* t = parseRawTreeNode(tree_string, root);
	return(t);
}

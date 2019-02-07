#include "TreeParser.h"

#include "../Environment.h"
#include "../IO.h"

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
    //std::cout << "Char: " << character << " Depth: " << tree_depth << std::endl;
    if (character == '(')
      tree_depth++;
    else if (character == ')')
      tree_depth--;
    else if (character == ',' and tree_depth == 0) {
      b.first = branch_string.substr(0, position);
      b.second = branch_string.substr(position + 1, branch_string.size() - position + 1);
      return(b);
    }
  }
  // Node with a single decendent.
  b.first = branch_string;
  b.second = "";
  return(b);
}

node_data IO::deconstructNodeString(std::string node_string) {
  static int ID = 0;
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
    branch_strings = splitBranchString(node_string.substr(1, last_parens_position - 1));
  }
  node_data n = {name, distance, branch_strings.first, branch_strings.second};
  return(n);
}

IO::RawTreeNode* IO::parseRawTreeNode(std::string node_string, RawTreeNode* up) {
  node_data n = deconstructNodeString(node_string);

  std::cout << node_string << std::endl;
  IO::RawTreeNode* t = new RawTreeNode;
  IO::RawTreeNode* left;
  IO::RawTreeNode* right;

  left = (n.left != "") ? parseRawTreeNode(n.left, t) : 0;
  right = (n.right != "") ? parseRawTreeNode(n.right, t) : 0;

  std::cout << "Name: " << n.name << " Distance: " << n.distance << std::endl;
  *t = {n.name, n.distance, up, left, right};
  return(t);
}

IO::RawTreeNode* IO::parseTree(std::string tree_string) {
	tree_string = cleanTreeString(tree_string);
	IO::RawTreeNode* root = new RawTreeNode;
	*root = {"Root", 0, NULL, NULL, NULL};
	IO::RawTreeNode* t = parseRawTreeNode(tree_string, root);
	return(t);
}

std::list<std::string> IO::getRawTreeNodeNames(IO::RawTreeNode* node) {
  std::list<std::string> node_names;
  node_names.push_back(node->name);
  if(node->left != 0) {
    node_names.splice(node_names.begin(), getRawTreeNodeNames(node->left));
  }
  if(node->right != 0) {
    node_names.splice(node_names.begin(), getRawTreeNodeNames(node->right));
  }
  return(node_names);
}

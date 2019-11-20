#include "TreeParser.h"

#include "../Environment.h"
#include "../IO/Files.h"

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

std::pair<std::string, std::string> IO::separate_node_info(std::string node_string) {
  /*
   * Takes a node string of form (...)name:distance or name:distance and returns a pair of strings.
   * The first element represents the potential child nodes (...), and the second element represents
   * the info name:distance.
   */
  int last_parens_position = node_string.find_last_of(')');
  if (last_parens_position == (signed int)string::npos) {
    std::cout << "Child nodes: " << "" << " Info: " << node_string << std::endl;
    return(std::pair<std::string, std::string>("", node_string));
  } else {
    std::cout << "Child nodes: " << node_string.substr(0, last_parens_position + 1) << " Info: " << node_string.substr(last_parens_position + 1) << std::endl;
    return(std::pair<std::string, std::string>(node_string.substr(0, last_parens_position + 1),
					     node_string.substr(last_parens_position + 1)));
  }
}

node_info IO::parse_node_info(std::string info_string) {
  static int ID = 0;
  int last_colon_position = info_string.find_last_of(':');
  std::string name;
  float distance = atof(info_string.substr(last_colon_position + 1).c_str());
  float bootstrap =  0.0;
  if(last_colon_position == 0) {
    name = "BNode" + std::to_string(ID);
    ID++;
  } else {
    std::string name_or_bootstrap = info_string.substr(0, last_colon_position);
    if(std::atof(name_or_bootstrap.c_str()) != 0.0) {
      bootstrap = std::atof(name_or_bootstrap.c_str());
      name = "BNode" + std::to_string(ID);
      ID++;
    } else {
      name = name_or_bootstrap;
    }
  }
  return(node_info({name, distance, bootstrap}));
}

std::vector<std::string> IO::separate_node_string(std::string node_string) {
  int depth = 0;
  std::vector<std::string>sub_nodes = {};

  if(node_string == "") {
    return(sub_nodes);
  }

  unsigned int previous_comma = 0;
  for(unsigned int i = 0; i < node_string.size(); i++) {
    char c = node_string[i];
    switch(c) {
    case '(':
      depth++;
      break;
    case ')':
      depth--;
      break;
    case ',':
      if(depth == 1) {
	sub_nodes.push_back(node_string.substr(previous_comma + 1, i - previous_comma - 1));
	previous_comma = i;
      }
    }
  }

  sub_nodes.push_back(node_string.substr(previous_comma + 1, node_string.size() - previous_comma - 2));

  //std::cout << node_string << std::endl;
  //for(auto it = sub_nodes.begin(); it != sub_nodes.end(); ++it) {
  //  std::cout << "-> " << *it << std::endl;
  // }
  //std::cout << std::endl;
  return(sub_nodes);
}

IO::RawTreeNode* IO::parseRawTreeNode(std::string node_string, RawTreeNode* up) {
  std::pair<std::string, std::string> node_pair = separate_node_info(node_string);
  node_info info = parse_node_info(node_pair.second);
  std::vector<std::string> child_nodes = separate_node_string(node_pair.first);
  std::cout << "Name: " << info.name << " Dist: " << info.distance << " Bootstrap: " << info.bootstrap << std::endl;

  IO::RawTreeNode* t = new RawTreeNode;
  IO::RawTreeNode* left;
  IO::RawTreeNode* right;

  switch(child_nodes.size()) {
  case 0:
    left = 0;
    right = 0;
    break;
  case 1:
    std::cerr << "Error: I have not implimented a tree parser to deal with internal nodes yet." << std::endl;
    exit(EXIT_FAILURE);
    left = parseRawTreeNode(child_nodes[0], t);
    right = nullptr;
    break;
  case 2:
    left = parseRawTreeNode(child_nodes[0], t);
    right = parseRawTreeNode(child_nodes[1], t);
    break;
  case 3:
    if(env.get<bool>("TREE.resolve_root")) {
      left = new RawTreeNode;
      *left = {"RRNode", 0.0000000001, t, parseRawTreeNode(child_nodes[0], left), parseRawTreeNode(child_nodes[1], left)};
      right = parseRawTreeNode(child_nodes[2], t);
    } else {
      std::cerr << "Error: detected an unrooted tree during tree parsing. Set TREE.resolve_root to arbitrarily root." << std::endl;
      exit(EXIT_FAILURE);
    }
    break;
  default:
    std::cerr << "Error: invalid number of downstream nodes during tree parsing." << std::endl;
    exit(EXIT_FAILURE);
    break;
  }

  *t = {info.name, info.distance * env.get<double>("TREE.scale_factor"), up, left, right};
  return(t);
}

IO::RawTreeNode* IO::parseTree(std::string tree_string) {
  tree_string = cleanTreeString(tree_string);
  std::cout << "Start: " << tree_string << std::endl;
  IO::RawTreeNode* root = new RawTreeNode;
  *root = {"Root", 0, NULL, NULL, NULL};
  IO::RawTreeNode* t = parseRawTreeNode(tree_string, root);

  if(t->distance > 0.0) {
    std::cout << "Warning: truncating the root branch node." << std::endl;
    t->distance = 0.0;
  };

  IO::printRawTree(t);
  //exit(EXIT_FAILURE);
  return(t);
}

// Utils

std::list<std::string> IO::getRawTreeNodeNames(const IO::RawTreeNode* node) {
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

void IO::printRawTree(const RawTreeNode* node) {
  std::cout << "Node: " << node->name << " " << node->distance << std::endl;
  if(node->left != nullptr) {
    printRawTree(node->left);
  }
  if(node->right != nullptr) {
    printRawTree(node->right);
  }
}

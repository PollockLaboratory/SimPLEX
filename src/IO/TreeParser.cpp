#include "TreeParser.h"

#include "../Environment.h"
#include "../IO/Files.h"

extern Environment env;
extern IO::Files files;

void print_node_info(const node_info& info) {
  std::cout << "[ " << info.name << " " << info.distance << " "
	    << info.bootstrap << " ]" << std::endl;
}

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
    return(std::pair<std::string, std::string>("", node_string));
  } else {
    return(std::pair<std::string, std::string>(node_string.substr(0, last_parens_position + 1),
					     node_string.substr(last_parens_position + 1)));
  }
}

bool is_number(const std::string& s) {
  auto it = s.begin();
  while(it != s.end() && std::isdigit(*it)) {
    ++it;
  }
  return(!s.empty() && it == s.end());
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
    if(is_number(name_or_bootstrap)) {
      bootstrap = std::atof(name_or_bootstrap.c_str());
      name = "BNode" + std::to_string(ID);
      ID++;
    } else {
      name = name_or_bootstrap;
    }
  }
  return(node_info({name, distance, bootstrap}));
}

std::queue<std::string> IO::separate_node_string(std::string node_string) {
  int depth = 0;
  std::queue<std::string>sub_nodes = {};

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
	sub_nodes.push(node_string.substr(previous_comma + 1, i - previous_comma - 1));
	previous_comma = i;
      }
    }
  }

  sub_nodes.push(node_string.substr(previous_comma + 1, node_string.size() - previous_comma - 2));

  return(sub_nodes);
}

std::pair<IO::RawTreeNode*, IO::RawTreeNode*> IO::resolve_child_nodes(std::queue<std::string> nodes, IO::RawTreeNode* up) {
  static int ID = 0;
  IO::RawTreeNode* left = parseRawTreeNode(nodes.front(), up);
  nodes.pop();
  if(nodes.size() == 1) {
    return(std::pair<IO::RawTreeNode*, IO::RawTreeNode*>(left, parseRawTreeNode(nodes.front(), up)));
  } else {
    if(env.get<bool>("TREE.resolve_multifurcation")) {
      RawTreeNode* intermediate = new RawTreeNode;
      std::pair<IO::RawTreeNode*, IO::RawTreeNode*> child_nodes = resolve_child_nodes(nodes, intermediate);
      *intermediate = {"MNode" + std::to_string(ID), 0.00000000001, up, child_nodes.first, child_nodes.second};
      ID++;
      return(std::pair<IO::RawTreeNode*, IO::RawTreeNode*>(left, intermediate));
    } else {
      std::cerr << "Error: detected nodes with multiple decendents. Set TREE.resolve_multifurcation to force resolution." << std::endl;
      exit(EXIT_FAILURE);
    }
  }
}

IO::RawTreeNode* IO::parseRawTreeNode(std::string node_string, RawTreeNode* up) {
  std::pair<std::string, std::string> node_pair = separate_node_info(node_string);
  node_info info = parse_node_info(node_pair.second);
  std::queue<std::string> child_nodes = separate_node_string(node_pair.first);
 
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
    left = parseRawTreeNode(child_nodes.front(), t);
    right = nullptr;
    break;
  default:
    std::pair<IO::RawTreeNode*, IO::RawTreeNode*> children = resolve_child_nodes(child_nodes, up);
    left = children.first;
    right = children.second;
    break;
  }
  *t = {info.name, info.distance, up, left, right};
  return(t);
}

IO::RawTreeNode* IO::parseTree(std::string tree_string) {
  tree_string = cleanTreeString(tree_string);
  IO::RawTreeNode* root = new RawTreeNode;
  *root = {"Root", 0, NULL, NULL, NULL};
  IO::RawTreeNode* t = parseRawTreeNode(tree_string, root);

  if(t->distance > 0.0) {
    std::cout << "Warning: truncating the root branch node." << std::endl;
    t->distance = 0.0;
  };
  
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

std::list<std::string> IO::getRawTreeNodeTipNames(const IO::RawTreeNode* node) {
  std::list<std::string> node_names;
  if(node->left == nullptr and node->right == nullptr) {
    node_names.push_back(node->name);
  }
  if(node->left != 0) {
    node_names.splice(node_names.begin(), getRawTreeNodeTipNames(node->left));
  }
  if(node->right != 0) {
    node_names.splice(node_names.begin(), getRawTreeNodeTipNames(node->right));
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

float IO::findRawTreeTotalLength(const RawTreeNode* node) {
  float total = node->distance;
  if(node->left != nullptr) {
    total += findRawTreeTotalLength(node->left);
  }
  if(node->right != nullptr) {
    total += findRawTreeTotalLength(node->right);
  }

  return(total);
}

void IO::deleteTree(RawTreeNode* node) {
  if(node->left != nullptr) {
    deleteTree(node->left);
  }

  if(node->right != nullptr) {
    deleteTree(node->right);
  }
}

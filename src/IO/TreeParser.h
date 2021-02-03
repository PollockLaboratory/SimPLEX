#ifndef TreeParser_h_
#define TreeParser_h_

#include <string> 
#include <iostream>
#include <list>
#include <utility>
#include <queue>

struct node_info {
  std::string name;
  float distance;
  float bootstrap; // The bootstrap value.
};

void print_node_info(const node_info&);

namespace IO {
  struct RawTreeNode {
    std::string name;
    double distance;
    RawTreeNode* up;
    RawTreeNode* left;
    RawTreeNode* right;
  };

  // Parsing the string/tokens of the newick file. 
  inline std::string& cleanTreeString(std::string &tree_string);
  std::pair<std::string, std::string> separate_node_info(std::string node_string);
  node_info parse_node_info(std::string info_string);
  std::queue<std::string> separate_node_string(std::string);

  // Helper functions for building the tree data structure.
  std::pair<RawTreeNode*, RawTreeNode*> resolve_child_nodes(std::queue<std::string>, RawTreeNode* up);
  RawTreeNode* parseRawTreeNode(std::string node_string, RawTreeNode* up);

  // Complete parse function.
  RawTreeNode* parseTree(std::string tree_string);

  // Utility functions.
  std::list<std::string> getRawTreeNodeNames(const RawTreeNode* node);
  std::list<std::string> getRawTreeNodeTipNames(const RawTreeNode* node);
  void printRawTree(const RawTreeNode* node);
  float findRawTreeTotalLength(const RawTreeNode* node);

  void deleteTree(RawTreeNode* node);
}

#endif



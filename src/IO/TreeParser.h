#ifndef TreeParser_h_
#define TreeParser_h_

#include <string> 
#include <iostream>
#include <list>
#include <vector>
#include <utility>

struct node_info {
  std::string name;
  float distance;
  float bootstrap; // The bootstrap value.
};

namespace IO {
  struct RawTreeNode {
    std::string name;
    double distance;
    RawTreeNode* up;
    RawTreeNode* left;
    RawTreeNode* right;
  };

  inline std::string& cleanTreeString(std::string &tree_string);
  std::pair<std::string, std::string> separate_node_info(std::string node_string);
  node_info parse_node_info(std::string info_string);
  std::vector<std::string> separate_node_string(std::string);

  RawTreeNode* parseRawTreeNode(std::string node_string, RawTreeNode* up);
  RawTreeNode* parseTree(std::string tree_string);
  std::list<std::string> getRawTreeNodeNames(const RawTreeNode* node);

  void printRawTree(const RawTreeNode* node);
}

#endif



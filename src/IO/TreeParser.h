#ifndef TreeParser_h_
#define TreeParser_h_

#include <string> 
#include <iostream>
#include <list>
#include <utility>

struct node_data {
  std::string name;
  float distance;
  std::string left;
  std::string right;
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
  node_data deconstructNodeString(std::string node_string);
  std::pair<std::string, std::string> splitBranchString(std::string branch_string);

  RawTreeNode* parseRawTreeNode(std::string node_string, RawTreeNode* up);
  RawTreeNode* parseTree(std::string tree_string);
  std::list<std::string> getRawTreeNodeNames(const RawTreeNode* node);

  void printRawTree(const RawTreeNode* node);
}

#endif



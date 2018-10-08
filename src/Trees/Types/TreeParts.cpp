/*
 * TREE PARTS
 * This file holds the classes for the components of the tree:
 * 	- Tree nodes - Contains the sequences.
 * 	- Branche segments - Contains the branch lengths and location of substitutions.
 */

#include "TreeParts.h"
#include <iostream>

// Branch segment.
BranchSegment::BranchSegment(float distance) {
	this->distance = distance;
//	std::cout << "Making new branch segment. Distance: " << this->distance << std::endl;
}

// Tree nodes.
TreeNode::TreeNode() {
	this->name = "NULL";
}

TreeNode::TreeNode(IO::RawTreeNode* raw_tree) {
	/* 
	 * TreeNode Constructor - from rawTreeNode.
	 */
	name = raw_tree->name;
	distance = raw_tree->distance;
}

TreeNode::TreeNode(std::string n) {
	/* 
	 * TreeNode Constructor - from name and distance.
	 */
	name = n;
}

bool TreeNode::isTip() {
	if(this->left == 0 and this->right == 0) {
		return(true);
	} else {
		return(false);
	}
}



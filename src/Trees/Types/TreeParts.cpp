/*
 * TREE PARTS
 * This file holds the classes for the components of the tree:
 * 	- Tree nodes - Contains the sequences.
 * 	- Branche segments - Contains the branch lengths and location of substitutions.
 */

#include "TreeParts.h"
#include <iostream>

extern double Random();

// Branch segment.
BranchSegment::BranchSegment(float distance) {
	this->distance = distance;
//	std::cout << "Making new branch segment. Distance: " << this->distance << std::endl;
}

BranchSegment::~BranchSegment() {
	decendant->up = 0;
}

// Tree nodes.

int TreeNode::unique_id = 0;

TreeNode::TreeNode(IO::RawTreeNode* raw_tree) {
	/* 
	 * TreeNode Constructor - from rawTreeNode.
	 */
	name = raw_tree->name;
	distance = raw_tree->distance;
	up = 0;
	left = 0;
	right = 0;

}

TreeNode::TreeNode(std::string n) {
	/* 
	 * TreeNode Constructor.
	 */
	name = n;
	up = 0;
	left = 0;
	right = 0;
}

TreeNode::TreeNode() {
	/* 
	 * TreeNode Constructor - from name and distance.
	 */
	name = "Node" + std::to_string(unique_id);
	unique_id++;
	up = 0;
	left = 0;
	right = 0;
}

bool TreeNode::isTip() {
	if(this->left == 0 and this->right == 0) {
		return(true);
	} else {
		return(false);
	}
}



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
	rates = {};
//	std::cout << "Making new branch segment. Distance: " << this->distance << std::endl;
}

BranchSegment::~BranchSegment() {
	decendant->up = 0;
}

std::ostream& operator<< (std::ostream &out, const BranchSegment &b) {
	out << b.distance;
	return out;
}

std::pair<int, int> BranchSegment::countSubstitutions() {
	/*
	 * The first int is number of 0 subs (C^0_b') and the second is the number of 1 sub (C^1_b')
	 */
	int num0subs = 0;
	int num1subs = 0;

	std::vector<int> s1 = *(ancestral->sequence);
	std::vector<int> s2 = *(decendant->sequence);
	for(int i = 0; i < s1.size(); i++) {
		if(s1[i] == s2[i]) {
			num0subs += 1;
		} else {
			num1subs += 1;
		}
	}
	
	std::pair<int, int> c = {num0subs, num1subs};
	return(c);
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



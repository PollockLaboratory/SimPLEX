#ifndef TreeParts_h_
#define TreeParts_h_

#include <string>
#include <map>
#include <vector>

#include "Trees/TreeParser.h"

using std::string;
using std::map;
using std::vector;

class TreeNode;

class BranchSegment {
	public:
		float distance;
		TreeNode* ancestral;
		TreeNode* decendant;

		BranchSegment(float distance);
};

class TreeNode {
	public:
		std::string name;
		double distance;
		vector<int> sequence;	
		BranchSegment* up;
		BranchSegment* left;
		BranchSegment* right;
		
		TreeNode();
		TreeNode(IO::RawTreeNode* raw_tree);
		TreeNode(std::string n);

		bool isTip();
};

#endif

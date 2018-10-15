#ifndef TreeParts_h_
#define TreeParts_h_

#include <string>
#include <map>
#include <vector>

#include "Trees/TreeParser.h"
#include "Sequence.h"

using std::string;
using std::map;
using std::vector;

class TreeNode;

class BranchSegment {
	public:
		float distance;
		TreeNode* ancestral;
		TreeNode* decendant;
		std::list<substitution> subs;

		BranchSegment(float distance);
		~BranchSegment();
};

class TreeNode {
	public:
		static int unique_id;
		std::string name;
		double distance;
		Sequence* sequence;	
		BranchSegment* up;
		BranchSegment* left;
		BranchSegment* right;
		
		TreeNode();
		TreeNode(IO::RawTreeNode* raw_tree);
		TreeNode(std::string n);

		bool isTip();
};
#endif

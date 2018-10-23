#ifndef TreeParts_h_
#define TreeParts_h_

#include <string>
#include <map>
#include <vector>
#include <iostream>

#include "Trees/TreeParser.h"
#include "Sequence.h"
#include "RateVector.h"

using std::string;
using std::map;
using std::vector;

class TreeNode;

class BranchSegment {
	public:
		BranchSegment(float distance);
		~BranchSegment();

		float distance;
		TreeNode* ancestral;
		TreeNode* decendant;
		std::list<substitution> subs;
		std::vector<RateVector*> rates;

		friend std::ostream& operator<< (std::ostream &out, const BranchSegment &b);

		std::pair<int, int> countSubstitutions();
};

class TreeNode {
	public:
		static int unique_id;
		std::string name;
		double distance;
		std::vector<int>* sequence;	
		BranchSegment* up;
		BranchSegment* left;
		BranchSegment* right;
		
		TreeNode();
		TreeNode(IO::RawTreeNode* raw_tree);
		TreeNode(std::string n);

		bool isTip();
};
#endif

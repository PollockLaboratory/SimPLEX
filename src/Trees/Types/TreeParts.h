#ifndef TreeParts_h_
#define TreeParts_h_

#include <string>
#include <map>
#include <vector>
#include <iostream>
#include <math.h>

#include "Trees/TreeParser.h"
#include "Sequence.h"
#include "RateVector.h"
#include "SubstitutionModels/SubstitutionModel.h"

class TreeNode;

class BranchSegment {
	public:
		BranchSegment(float distance);
		~BranchSegment();

		float distance;
		TreeNode* ancestral;
		TreeNode* decendant;

		void set_rate_vector(int pos);
		void set_rate_vector(int pos, RateVector* rv);
		double get_rate(int pos, int dec_state); 			

		// Key statistics.
		std::vector<substitution> subs;
		int num0subs;
		int num1subs;

		bool virtualSubstituionQ(int state);
		void updateStats();

		friend std::ostream& operator<< (std::ostream &out, const BranchSegment &b);
	private:
		std::vector<RateVector*> rates; //By site.
};

class TreeNode {
	public:
		static int unique_id;
		std::string name;
		double distance;
		
		BranchSegment* up;
		BranchSegment* left;
		BranchSegment* right;
		
		std::vector<int>* sequence;	
		SequenceAlignment* MSA;
	
		SubstitutionModel* SM;

		TreeNode();
		TreeNode(IO::RawTreeNode* raw_tree);
		TreeNode(std::string n);

		void sampleSequence();
		void sampleSinglePosition(int pos);

		bool isTip();
		bool sampled;
};

#endif

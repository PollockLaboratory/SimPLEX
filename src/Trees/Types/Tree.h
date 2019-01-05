#ifndef Tree_h_
#define Tree_h_

#include <algorithm>
#include <functional>
#include <istream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <list>

#include "Sequence.h"
#include "SubstitutionModel.h"
#include "Trees/TreeParser.h"
#include "TreeParts.h"
#include "BranchSplitting.h"

using std::string;
using std::map;
using std::vector;

class Tree {
	public:
		TreeNode* root;
		SequenceAlignment* MSA;
		SubstitutionModel* SM;
		std::map<float, std::pair<int, int>> substitution_counts;

		Tree();
		Tree& operator=(Tree tree);

		int seqLen;
		map<string, vector<int>> names_to_sequences;
		std::list<BranchSegment*> branchList; // Potentially should be vectors.
		std::vector<TreeNode*> nodeList;

		// Internal Tree nodes.
		void connectNodes(TreeNode* &ancestral, BranchSegment* &ancestralBP, TreeNode* &decendant, float distance);
		TreeNode* createTreeNode(IO::RawTreeNode* raw_tree, TreeNode* &ancestralNode, BranchSegment* &ancestralBP);

		// Setting options.
		float max_seg_len; // Max segment length.
		std::function< std::pair<BranchSegment*, BranchSegment*>(float)> splitBranch; // Algorithm for splitting branches.
		
		// Initializing.
		void Initialize(IO::RawTreeNode* raw_tree, SequenceAlignment* &MSA, SubstitutionModel* &SM);

		// Debug tools.
		void printBranchList();
		void printNodeList();
		void printParameters();

		// Sampling.
		bool SampleParameters();
		
		// Recording state data.
		void RecordTree();
		void RecordState(int gen, double l);
	
		// Likelihood.
		void findKeyStatistics(); //Find the key statistics need for the likelihood function.
		double calculate_likelihood();
		double partial_calculate_likelihood();
	private:
		void configureSequences(TreeNode* n);
		void configureRateVectors();
		
		float u;
	public:
//		typedef unsigned int size_type;
		int id;

		static std::ofstream tree_out;
		static std::ofstream substitutions_out;

		void InitializeOutputStreams();
};

#endif

#ifndef Tree_B1_h_
#define Tree_B1_h_

#include "Tree.h"

using std::string;
using std::map;
using std::vector;

class Tree_B1 : public Tree {

//	friend class Model;
public:

	/* Inherited properties
	string name;
	double distance;
	vector<int> sequence;
	*/
	
	
	//Tree_B1* left;
	//Tree_B1* right;
	//Tree_B1* up;
	

	Tree_B1();
	Tree_B1(const Tree_B1& tree);
	Tree_B1& operator=(Tree_B1 tree);
	// ~Tree_B1();

	Tree_B1* Clone();
	
	void Initialize(map<string, vector<int> > taxa_names_to_sequences,
			vector<string> states);

	void SampleParameters();

	//void RecordState();

	//bool IsSubtree() const;
	bool IsSegment() const;
	//bool IsLeaf() const;


protected:
	/* Inherited
	typedef unsigned int size_type;

	static int number_of_trees;
	int id;
	bool is_constant;



	vector<string> states;



	static std::ofstream tree_out;
	static std::ofstream substitutions_out;
	static std::ofstream sequences_out;
	*/
	
	//void ReadFromTreeFile();
	void ReadFromString(string);
	void ReadSubtree(string);
	//void ExtractDistance(string&);
	void SegmentBranch();


	void InitializeSequences(
			std::map<string, vector<int> > taxa_names_to_sequences);
	//void InitializeOutputStreams();

	void SampleSubtreeParameters();
	//void SampleSequence();
	//void SampleDistance();

	void RecordSubtreeState();
	//void RecordSequence();
	//void RecordToTreefile();
	void RecordSubstitutions();
	//void RecordChildSubstitutions(Tree* child);
	//void AddGenerationEndIndicatorsToOutputFiles();

	//string IdToString();
	//void DescendentStateSampling();
	//void MetropolisHastingsStateSampling();
};

#endif

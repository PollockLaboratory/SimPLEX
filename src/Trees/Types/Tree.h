#ifndef Tree_h_
#define Tree_h_

#include <algorithm>
#include <istream>
#include <fstream>
#include <string>
#include <map>
#include <vector>

using std::string;
using std::map;
using std::vector;

class Tree {

//	friend class Model;
public:

	string name;
	double distance;
	vector<int> sequence;
	// STP: How can I allow a derived class to put a pointer to the derived
	// class in these members?
	Tree* left;
	Tree* right;
	Tree* up;

	Tree();
	Tree(const Tree& tree);
	Tree& operator=(Tree tree);
	virtual ~Tree();
	
	//STP: This method is necessary for recursive inherited classes to work
	virtual Tree* Clone();
	virtual void Initialize(map<string, vector<int> > taxa_names_to_sequences,
			vector<string> states);
	virtual void SampleParameters();
	virtual void RecordState();

	bool IsRoot() const;
	virtual bool IsSubtree() const;
	virtual bool IsLeaf() const;

// STP: All these things should be protected, but when they are, the derived 
// class Tree_B1 seems to not have access to them. The definition of protected
// here means derived classes should have access to them. 
//protected:
public:
	typedef unsigned int size_type;

	static int num_trees;
	int id;
	bool is_constant;


	/* Should this really be _part of_ the tree? Or should the tree simply know
	 * about it?
	 * Maybe this should be a member of the model. Or perhaps this should
	 * simply be a pointer to the integer_to_state.
	 *
	 * No. The tree should have everything it needs to print itself within
	 * itself. It should have all the information in it. That means integer_
	 * to_state must be a member of every tree even though it is identical
	 * along all trees.... hm maybe it should be a class static then.
	 * Or I could use a shared pointer! That way there is only one matrix. And
	 * different instances of trees can have different outputs.
	 *
	 * This also guarantees that the object will always exist if there is a
	 * shared pointer to it (unlike a regular pointer).
	 */
 
	vector<string> states;
	
	
	/**
	 * Should these really be class statics? I like pointers better than class
	 * statics. And now since there will be a terminate call, it can delete
	 * anything allocated on the heap. So perhaps a better solution would be
	 * to allocate the ofstreams on the heap at initialization and then make
	 * pointers to them and delete them at termination.
	 *
	 * Or use shared pointers instead of manually allocating them.
	 *
	 */


	static std::ofstream tree_out;
	static std::ofstream substitutions_out;
	static std::ofstream sequences_out;

	void ReadFromTreeFile();
	virtual void ReadFromString(string);
	virtual void ReadSubtree(string);
	void ExtractDistance(string&);
	void ExtractName(string&);


	virtual void InitializeSequences(
			std::map<string, vector<int> > taxa_names_to_sequences);
	void InitializeOutputStreams();

	virtual void SampleSubtreeParameters();
	void SampleSequence();
	void SampleDistance();

	virtual void RecordSubtreeState();
	void RecordSequence();
	void RecordToTreefile();
	virtual void RecordSubstitutions();
	void RecordChildSubstitutions(Tree* child);
	void AddGenerationEndIndicatorsToOutputFiles();

	string IdToString();
	virtual void DescendentStateSampling();
	void MetropolisHastingsStateSampling();
};

#endif

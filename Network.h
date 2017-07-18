#ifndef Network_h_
#define Network_h_

#include <algorithm>
#include <istream>
#include <fstream>
#include <string>
#include <map>
#include <vector>

#include "Matrix.h"

using std::string;
using std::map;
using std::vector;


struct Node {
	string name;
	vector<int> sequence;
};

class Network {

	friend class Model;
public:
	Network();
	Network(const Network& network);

	void Initialize(map<string, vector<int> > taxa_names_to_sequences, map<int, string> integer_to_state);

	void SampleParameters();

	void RecordState();

private:
	typedef unsigned int size_type;

	static int number_of_networks;
	int id;
	bool is_constant;

	vector<Node> nodes;
	Matrix<double> distances;


	//Should this really be _part of_ the Network? Or should the Network simply know
	// about it?
	// Maybe this should be a member of the model. Or perhaps this should
	// simply be a pointer to the integer_to_state
	map<int, string> integer_to_state;

	/**
	 * Should these really be class statics? I like pointers better than class
	 * statics. And now since there will be a terminate call, it can delete
	 * anything allocated on the heap. So perhaps a better solution would be
	 * to allocate the ofstreams on the heap at initialization and then make
	 * pointers to them and delete them at termination.
	 *
	 */
	static std::ofstream Network_out;
	static std::ofstream substitutions_out;
	static std::ofstream sequences_out;

	void InitializeSequences(std::map<string, vector<int> > taxa_names_to_sequences);
	void InitializeOutputStreams();

	void SampleSequences();

	void RecordSubNetworkState();
	void RecordSequence();
	void RecordNameAndDistanceToNetworkfile();
	void RecordSubstitutions();
	void RecordChildSubstitutions(Network* child);
	void AddGenerationEndIndicatorsToOutputFiles();
	void SampleSubNetworkParameters();
	std::string IdToString();
	void SampleDistances();
};

#endif

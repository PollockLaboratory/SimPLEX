#ifndef Sequence_h_
#define Sequence_h_

#include<string>
#include<map>
#include<vector>
#include<list>
#include<set>

struct substitution {
	int pos;
	int anc;
	int dec;
};

class Sequence;

class SequenceAlignment {
	public:
		SequenceAlignment();
		std::map<std::string, Sequence*> taxa_names_to_sequences;
		
		std::set<int> columns_with_gaps;
		std::vector<int> columns_without_gaps;

		std::vector<std::string> states;
		std::map<std::string, int> state_to_integer;
		std::map<int, std::string> integer_to_state;

		// Adding sequences to alignment.
		void add(std::string name, std::string sequence_str);
		void addVariable(std::string name);

		// Processing input sequences.	
		std::vector<int> EncodeSequence(std::string sequence);
		void AddStateToStates(std::string state);
		void DetermineColumnsWithoutGaps();
		void RemoveColumnsWithGapsFromSequences();
		std::vector<int> RemoveGapsFromEncodedSequence(std::vector<int> encoded_sequence);

		void Initialize();
		void print();

		// Utilities
		static std::vector<int> findParsimony(const std::vector<int> &s1, const std::vector<int> &s2);
		static std::list<substitution> findSubstitutions(const std::vector<int> &anc, const std::vector<int> &dec);
};

class Sequence {
	public: 

		std::vector<int> encoded_sequence;
		SequenceAlignment* MSA;
		
		Sequence(std::vector<int> enc, SequenceAlignment* MSA);
		void set(std::vector<int> seq);
		std::string as_str();
};

#endif

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
		std::map<std::string, std::vector<int>> taxa_names_to_sequences;
		
		std::set<int> columns_with_gaps;
		std::vector<int> columns_without_gaps;

		std::vector<std::string> states;
		std::map<std::string, int> state_to_integer;
		std::map<int, std::string> integer_to_state;

		// Adding sequences to alignment.
		void add(std::string name, std::string sequence_str);
		void add(std::string name);

		// Processing input sequences.	
		std::vector<int> EncodeSequence(std::string sequence);
		void AddStateToStates(std::string state);
		void DetermineColumnsWithoutGaps();
		void RemoveColumnsWithGapsFromSequences();
		std::vector<int> RemoveGapsFromEncodedSequence(std::vector<int> encoded_sequence);

		void Initialize();
		void print();

		// Utilities
		int numCols();
		std::string decodeChar(int &c);
		std::string decodeSequence(std::vector<int> &enc_seq);
		static std::vector<int> findParsimony(const std::vector<int> &s1, const std::vector<int> &s2);
		static std::list<substitution> findSubstitutions(const std::vector<int> &anc, const std::vector<int> &dec);
};
#endif

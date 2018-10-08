#ifndef Sequence_h_
#define Sequence_h_

#include<string>
#include<map>
#include<vector>
#include<set>

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

		void add(std::string name, std::string sequence_str);
		std::vector<int> EncodeSequence(std::string sequence);
		void AddStateToStates(std::string state);
		void DetermineColumnsWithoutGaps();
		void RemoveColumnsWithGapsFromSequences();
		std::vector<int> RemoveGapsFromEncodedSequence(std::vector<int> encoded_sequence);

		void Initialize();
		void print();
};

class Sequence {
	public: 

		std::vector<int> encoded_sequence;
		SequenceAlignment* MSA;
		
		Sequence(std::vector<int> enc, SequenceAlignment* MSA);
		std::string as_str();
};

#endif

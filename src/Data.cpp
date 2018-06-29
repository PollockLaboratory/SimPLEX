#include "Data.h"

#include <algorithm>
#include <map>

#include "Options.h"

extern Options options;

static const int gap_indicator = -1;

Data::Data() {
//	cout << "Data constructor" << endl;
}

Data::~Data() {
//	cout << "Data destructor" << endl;
}

void Data::Initialize() {
	ReadSequences();
	DetermineColumnsWithoutGaps();
	RemoveColumnsWithGapsFromSequences();
}

string Data::cleanLine(string line) {
	/* 
	 * Cleans up the a line to remove blank spaces a line returns.
	 */

	//Remove spaces and concatenate words
	line.erase(std::remove(line.begin(), line.end(), ' '), line.end());
	//This removes the carriage return in Windows files
	line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
	return line;
}

/// This function is too long.
void Data::ReadSequences() {
	// This is the only dependency on options.
	ifstream sequences_in(options.get("sequences_file").c_str());

	if (not sequences_in.good()) {
		cerr << "Cannot read sequence file" << endl;
		exit(-1);
	}

	vector<int> encoded_sequence;
	string line;
	string sequence = "";
	string name;

	while (sequences_in.good()) {
		getline(sequences_in, line);
		line = cleanLine(line);

		if (line == "")
			continue;

		if (line.at(0) == '>') {
			if (sequence != "") {
				taxa_names_to_sequences[name] = EncodeSequenceAndReportGaps(sequence);
			}
			name = line.substr(1);
			sequence = "";
		} else {
			sequence += line;
		}
	}
// Add final sequence
	taxa_names_to_sequences[name] = EncodeSequenceAndReportGaps(sequence);
}

vector<int> Data::EncodeSequenceAndReportGaps(string sequence) {
	vector<int> encoded_sequence(sequence.length());

	for (int site = 0; site < sequence.length(); site++) {
		string current_pos = sequence.substr(site, 1);
		// cout << state;
		if (HasGap(current_pos)) {
			ReportGapAtColumn(site);
			encoded_sequence.at(site) = gap_indicator;
		} else {
			AddStateToStates(current_pos);
			encoded_sequence.at(site) = state_to_integer[current_pos];
		}
	}
	// cout << endl;
	return encoded_sequence;
}

bool Data::HasGap(string sequence) {
	// When sequence.find == string::npos, the string was not found
	return (sequence.find('-') != string::npos
			or sequence.find('?') != string::npos);
}

void Data::ReportGapAtColumn(int column) {
	columns_with_gaps.insert(column);
}

void Data::AddStateToStates(string state) {
	// Checks to see of state exists in state_to_interger_map.
	if (state_to_integer.find(state) == state_to_integer.end()) {
		int encoded_state = state_to_integer.size();
		state_to_integer[state] = encoded_state;

		states.push_back(state);
	}
}

void Data::DetermineColumnsWithoutGaps() {
	int number_of_sites = taxa_names_to_sequences.begin()->second.size();
	cout << "Number of columns in alignment: " << number_of_sites << endl;

	for (int site = 0; site < number_of_sites; site++) {
		// If site is not in columns with gaps
		if (columns_with_gaps.find(site) == columns_with_gaps.end()) {
			columns_without_gaps.push_back(site);
		}
	}

	cout << "Number of columns without gaps: "
			<< columns_without_gaps.size() << endl;
}

void Data::RemoveColumnsWithGapsFromSequences() {
	if (options.debug)
		cout << "Removing gaps" << endl;

	for (std::map<string, vector<int> >::iterator it =
			taxa_names_to_sequences.begin();
			it != taxa_names_to_sequences.end(); it++) {
		vector<int> encoded_sequence = it->second;
		it->second = RemoveGapsFromEncodedSequence(encoded_sequence);
	}
}

vector<int> Data::RemoveGapsFromEncodedSequence(vector<int> encoded_sequence) {
	vector<int> encoded_sequence_without_gaps(columns_without_gaps.size());

	for (int site = 0; site < columns_without_gaps.size(); site++) {
		encoded_sequence_without_gaps.at(site) = encoded_sequence.at(
				columns_without_gaps.at(site));
	}
	return encoded_sequence_without_gaps;
}

void Data::PrintTaxaAndSequences() {
	for (std::map<string, vector<int> >::iterator it =
			taxa_names_to_sequences.begin();
			it != taxa_names_to_sequences.end(); ++it) {
		cout << it->first << ":";
		for (int i = 0; i < it->second.size(); i++)
			cout << " " << it->second.at(i);
		cout << endl;
	}
}

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include <algorithm>

#include "Matrix.h"

using namespace std;

void InitializeStateFromFile(std::string state_in_file, vector<string> states);

int main() {
	vector<string> states(4);
	states.at(0) = "A";
	states.at(1) = "C";
	states.at(2) = "G";
	states.at(3) = "T";

	InitializeStateFromFile("matrix", states);
}

void InitializeStateFromFile(std::string state_in_file, vector<string> states) {
	std::ifstream state_in(state_in_file.c_str());
	int number_of_pairs = states.size() * states.size();
	cout << "HI THERE" << endl;

// read header
	vector<string> labels(number_of_pairs);
	for (int pair = 0; pair < number_of_pairs; pair++) {

		state_in >> labels.at(pair);
		cout << labels.at(pair) << " ";
	}
	cout << endl;

	Matrix<double> matrix(states.size() , states.size());
	for (int label = 0; label < number_of_pairs; label++) {
		int row = find(states.begin(), states.end(),
				labels.at(label).substr(0, 1)) - states.begin();
		int column = find(states.begin(), states.end(), labels.at(label).substr(1,2))
						- states.begin();
		cout << row << " " << column << endl;
		state_in >> matrix.at(row, column);
	}

	for (int row = 0; row < matrix.number_of_rows; row++) {
		for (int column = 0; column < matrix.number_of_columns; column++) {
			cout << matrix.at(row, column) <<  " ";
		}
		cout << endl;
	}
}

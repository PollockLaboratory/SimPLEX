#ifndef Environment_h_
#define Environment_h_

#include <fstream>
#include <string>
#include <queue>
#include <map>

using std::string;

class Environment {
public:
	Environment();
	void PrintOptions();
	void ReadOptions(std::ifstream &default_file_stream, std::ifstream &options_file_stream);         // read in default and control files

	int seed;                   // random number generator seed
    	bool debug;                 // turns debugging on or off

	// Storing option values.
	int total_options;
	std::map<std::string, int> option_to_index;
	std::map<int, std::string> option_values;

	// Retreiving option values.
	void check_option_exists(std::string);
	std::string get(std::string);
	int get_int(std::string);
	float get_float(std::string);

private:
	void ReadControlFile(std::ifstream &file_stream);
	void SetOption(string option, string value);
	void InitializeRandomNumberGeneratorSeed();
};

#endif

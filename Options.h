

#ifndef Options_h_
#define Options_h_

#include <fstream>
#include <string>
#include <queue>


using std::string;

class Options {
public:
	Options();

	void ReadOptions();

	void PrependOutputDirectory(string &parameter);

	string default_control_file;
	string options_control_file;

	bool debug;

	int random_number_generator_seed;

	bool constant_tree;


	std::queue<int> substitution_model_types;
	std::queue<bool> constant_substitution_models;
	std::queue<int> substitution_models_initialization_type;
	std::queue<string> substitution_models_initialization_files;

	int number_of_substitution_models_in_mixture_model;


	int generations;
	int output_frequency;
	double max_segment_length;
	int tree_type;


	string output_directory;

	string tree_file;
	string sequences_file;
	string tree_out_file;
	string sequences_out_file;
	string substitutions_out_file;
	string likelihood_out_file;

private:
	void ReadControlFile(string controlfile);
	void ProcessOptions();
	void SetOption(string option, string value);
	//Might remove this function for lack of use
	void WriteOptions(std::ofstream& osparam);
	void UpdateOutputOptions();
	void CopyFile(string source_filename, string destination_filename);
	void InitializeRandomNumberGeneratorSeed();
};

#endif

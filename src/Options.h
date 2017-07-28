

#ifndef Options_h_
#define Options_h_

#include <fstream>
#include <string>
#include <queue>


using std::string;

class Options {
public:
	Options();
	void ReadOptions();         // read in default and control files
    	void PrependOutputDirectory(string &parameter);     // wtf?

	string defaultfile;         // where to find default settings
	string optionsfile;         // where to find optional control settings

	int seed;                   // random number generator seed
    	bool debug;                 // turns debugging on or off
	bool constant_tree;         // allow branch lengths, topology to vary?
    	int mixture_classes;        // number of classes in mixture model

	std::queue<int> substitution_model_types;
	std::queue<bool> constant_substitution_models;
	std::queue<int> substitution_models_initialization_type;
	std::queue<string> substitution_models_initialization_files;

	int gens;                   // mcmc generation end point
	int outfreq;                // how often in mcmc to print
	double max_segment_length;         // bigger than this and segments should be split
	int tree_type;              // we might get rid of this


	string outdir;              // directory for output files

	string treefile;            // name of file containing tree or trees
	string seqfile;             // name of file containing sequences
	string treeout;             // name of file to output trees
	string seqsout;             // name of file to output sequences
	string subsout;             // name of file to output substitutions
	string lnlout;              // name of file to output likelihoods

private:
	void ReadControlFile(string controlfile);
	void ProcessOptions();
	void SetOption(string option, string value);
	//Might remove this function for lack of use
	void WriteOptions(std::ofstream& osparam);
	void UpdateOutputOptions();
	void CopyFile(string source_filename, string destination_filename);
	void debugint(bool debug, string blurb, int integer);
	void InitializeRandomNumberGeneratorSeed();
};

#endif

#include "Options.h"

#include <cstdlib> // For atoi()
#include <iostream>

#ifdef _WIN32
#include <dir.h> // For mkdir()
#include <ctime> // For time()
//#include <Windows.h> // For CreateDirectory. Can't use because of CopyFile
#else
#include <sys/stat.h>// For making directories in Linux and OS X
#include <sys/times.h> // For time()
#endif

//Initialize all options here but don't do work.
Options::Options() {

	defaultfile = "resources/defaults.ctrl";   // default control file, sets the type of analysis
	optionsfile = "resources/options.ctrl";   // options control file, tend to vary from run to run
    	outdir = "testing/test_run/";              // directory for output files
    	treefile = "testing/5taxon.tree";           // name of file containing tree or trees
    	seqfile = "testing/5taxon.fasta";      // name of file containing sequences
    	treeout = "tree_out.newick";       // name of file to output trees
    	seqsout = "sequences_out.fasta";  // name of file to output sequences
    	subsout = "substitutions";  // name of file to output substitutions
    	lnlout = "likelihoods";  // name of file to output likelihoods

    
	seed = 830472018;               // random number generator seed
	debug = false;                  // factory default is don't debug
	constant_tree = false;          // factory default is that tree varies
	mixture_classes = 0;            // factory default is no mixtures
    	// model queues do not have default values

	gens =  100;                    // factory default is 100 generations
	outfreq = 25;                   // factory default is print output every 25 generations
	max_segment_length = 0.08;             // bigger than this and segments should be split
	tree_type = 0;                  // numbered tree types; should this be a name?
}

void Options::ReadOptions() {
	ReadControlFile(defaultfile);
	ReadControlFile(optionsfile);
	ProcessOptions();
}

void Options::ProcessOptions() {
	InitializeRandomNumberGeneratorSeed(); // InitRandom(seed);
	UpdateOutputOptions();                  // meaning what?

	this->CopyFile(optionsfile, outdir + optionsfile);
	this->CopyFile(defaultfile, outdir + defaultfile);

	if (substitution_model_types.size() > 1) { // Remove substitution model type read from default file
		substitution_model_types.pop();  // this seems clunky
	}  if (mixture_classes < 1) mixture_classes = 1;
}

void Options::ReadControlFile(string control_file) {
	std::cout << "Reading options from " << control_file << std::endl;

	std::ifstream controlfile_stream(control_file.c_str());
	if (not controlfile_stream.good()) {
		std::cerr << "Cannot read control file \"" << control_file << "\"" << std::endl;
		exit(-1);
	}

	bool in_comment = false;
	string key = "";  string value = "";
	while (controlfile_stream.good()) {   controlfile_stream >> key;
		if (key == "#") { in_comment = not in_comment;
		} else if (not in_comment) { // read value and set option
			controlfile_stream >> value;
			SetOption(key, value);
		}
	}
}

// there is no check if options are of right type
void Options::SetOption(string option, string value) {  // can we do a case here?
	if (option == "debug") debug = atoi(value.c_str());
	else if (option == "seed") seed = atoi(value.c_str());
	else if (option == "treefile") treefile = value;
	else if (option == "seqfile") seqfile = value;
	else if (option == "treeout") treeout = value;
	else if (option == "subsout") subsout = value;
	else if (option == "seqsout") seqsout = value;
	else if (option == "lnlout") lnlout = value;
	else if (option == "outfreq") outfreq = atoi(value.c_str());
	else if (option == "gens") gens = atoi(value.c_str());
	else if (option == "constant_tree") constant_tree = atoi(value.c_str());
	else if (option == "tree_type") tree_type = atoi(value.c_str());
	else if (option == "max_segment_length") max_segment_length = atof(value.c_str());
	else if (option == "substitution_model_type") substitution_model_types.push(atoi(value.c_str()));
	else if (option == "constant_substitution_model") constant_substitution_models.push(atoi(value.c_str()));
	else if (option == "substitution_model_initialization") substitution_models_initialization_type.push(atoi(value.c_str()));
	else if (option == "substitution_model_initialization_file") substitution_models_initialization_files.push(value);
	else if (option == "mixture_classes") mixture_classes = atoi(value.c_str());
	else if (option == "outdir") outdir = value;
	else { //STP: Unrecognized option is a non-fatal error so simply print warning
		std::cerr << "Unrecognized option found in control file:" << std::endl;
		std::cerr << "\"" << option << "\" with a value \"" << value << "\""
				<< std::endl;
	}
}

// Might remove this function because it is never called
// but it should be
void Options::WriteOptions(std::ofstream& out_stream) {
	/* Modified 2013-4-1 STP
	 * Using setw(20) cannot handle longer arguments well. In particular, the
	 * output files are long strings and are being cut at 20 characters.
	 * I changed setw to using \t instead
	 */
	out_stream << "Using options specified by control file:" << std::endl;
	out_stream << std::endl << "# debug and loudness controls #" << std::endl;
	out_stream << "\tdebug\t\t\t" << debug << std::endl;
	out_stream << "\tseed\t\t\t" << seed << std::endl;

	out_stream << std::endl << "# control generations, frequency of sampling #"
			<< std::endl;
	out_stream << "\tgenerations\t\t\t" << gens << std::endl;

	out_stream << std::endl
			<< "# options for simulation to get proposal step in MCMC #"
			<< std::endl;

	out_stream << std::endl << "# probablity calculation method #" << std::endl;

	out_stream << std::endl << "# optional input file names #" << std::endl;
	out_stream << "\ttreefile\t\t\t" << treefile << std::endl;
	out_stream << "\tseqfile\t\t\t" << seqfile << std::endl;

	out_stream << std::endl << "# optional output file names #" << std::endl;
	out_stream << "\ttreeoutfile\t\t\t" << treeout << std::endl; // Validated
	out_stream << "\tsuboutfile\t\t\t" << subsout << std::endl; // Validated
	out_stream << "\tseqoutfile\t\t\t" << seqsout << std::endl;
	out_stream << "\tlikelihoodoutfile\t\t\t" << lnlout
			<< std::endl;
}

void Options::UpdateOutputOptions() {
	if (outdir == "") {
		time_t rawtime = time(NULL); // initializes rawtime
		struct tm * timeinfo = localtime(&rawtime); // transforms it into a local tm struct

		char dir[100];
		strftime(dir, 100, "%Y-%m-%d_%H%M%S/", timeinfo);
		outdir = std::string(dir);
	}
	char lastchar = outdir.at(outdir.length() - 1);
	if (debug)
		std::cout << "The last char in " << outdir << " is "
				<< lastchar << std::endl;

	if (lastchar != '/' && lastchar != '\\') {
		if (debug)
			std::cout << "last char is not /" << std::endl;
		outdir += '/';
	}
	// For Windows
	// According to http://sourceforge.net/p/predef/wiki/OperatingSystems/
	// TODO: double check this works on windows for directory that already exists
#ifdef _WIN32
//	if (CreateDirectory(outputDirectory.c_str(), 0)) {
	if (mkdir(outdir.c_str())) {
		std::cout << "Could not make output directory " << outdir
				<< std::endl;
	} else {
		std::cout << "Making directory successful" << std::endl;
	}
#else	//ifdef __linux__ || __APPLE__
	struct stat st;
	if (stat(outdir.c_str(), &st) == 0 && S_ISDIR(st.st_mode)) {
		if (not debug) {
			std::cout << "output dir " << outdir
			<< " exists. Overwrite? (Y/n)" << std::endl;
			if (getchar() != 'Y')
			exit(1);
		}
	} else {
		mkdir(outdir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
	}
	/*Read, write, and search, or execute, for the file owner; S_IRWXU is the bitwise inclusive OR of S_IRUSR, S_IWUSR, and S_IXUSR.
	 Read, write, and search or execute permission for the file's group. S_IRWXG is the bitwise inclusive OR of S_IRGRP, S_IWGRP, and S_IXGRP.
	 Read, write, and search or execute permission for users other than the file owner. S_IRWXO is the bitwise inclusive OR of S_IROTH, S_IWOTH, and S_IXOTH.*/
//		std::cout << "making directory successful: " << m_strOutputDirectory << std::endl;
#endif

	PrependOutputDirectory(treeout);
	PrependOutputDirectory(subsout);
	PrependOutputDirectory(lnlout);
	PrependOutputDirectory(seqsout);
}

void Options::PrependOutputDirectory(std::string &parameter) {
	if (parameter == "") std::cerr << "Cannot prepend output directory to empty parameter" << std::endl;
	parameter = outdir + parameter;
}

void Options::debugint(bool debug, string blurb, int integer) {
    if (debug) std::cout << blurb << integer << std::endl;
}

void Options::InitializeRandomNumberGeneratorSeed() {
    string blurb = "Seed not specified so set to ";
    if (seed == 0) { // If seed not specified, set from clock
	    seed = time(0); 
	    debugint(debug,blurb,seed);
    }
    srand(seed);
}

void Options::CopyFile(string sourcefile, string newfile) {
	std::ifstream source(sourcefile.c_str());
	std::ofstream destination(newfile.c_str());
	destination << source.rdbuf();
}

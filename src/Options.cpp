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

	default_control_file = "options_default.ctrl";
	options_control_file = "options.ctrl";

	random_number_generator_seed = 0;

	debug = false;

	constant_tree = false;


	number_of_substitution_models_in_mixture_model = 0;



	generations = 0;

	output_frequency = 25;

	max_segment_length = 0.08;
	tree_type = 0;
}

void Options::ReadOptions() {
	ReadControlFile(default_control_file);
	ReadControlFile(options_control_file);
	ProcessOptions();
}

void Options::ProcessOptions() {
	InitializeRandomNumberGeneratorSeed();
	UpdateOutputOptions();

	this->CopyFile(options_control_file,
			output_directory + options_control_file);
	this->CopyFile(default_control_file,
			output_directory + default_control_file);

	if (substitution_model_types.size() > 1) {
		// Remove the substitution model type read from the default control
		// file
		substitution_model_types.pop();
	}

	if (number_of_substitution_models_in_mixture_model < 1)
		number_of_substitution_models_in_mixture_model = 1;
}

void Options::ReadControlFile(string control_file) {
	std::cout << "Reading options from " << control_file << std::endl;

	std::ifstream controlfile_stream(control_file.c_str());
	if (not controlfile_stream.good()) {
		std::cerr << "Cannot read control file \"" << control_file << "\""
				<< std::endl;
		exit(-1);
	}

	bool in_comment = false;
	string key = "";
	string value = "";

	while (controlfile_stream.good()) {
		controlfile_stream >> key;

		if (key == "#") {
			in_comment = not in_comment;
		} else if (not in_comment) {
			controlfile_stream >> value;
			SetOption(key, value);
		}
	}
}

void Options::SetOption(string option, string value) {
	if (option == "debug")
		debug = atoi(value.c_str());
	else if (option == "random_number_generator_seed") {
		random_number_generator_seed = atoi(value.c_str());
		if (debug)
			std::cout << "FOUND A SEED! The seed is "
					<< random_number_generator_seed << std::endl;
	} else if (option == "tree_file") {
		tree_file = value;
	} else if (option == "sequences_file")
		sequences_file = value;
	else if (option == "tree_out_file")
		tree_out_file = value;
	else if (option == "substitutions_out_file")
		substitutions_out_file = value;
	else if (option == "sequences_out_file")
		sequences_out_file = value;
	else if (option == "likelihood_out_file")
		likelihood_out_file = value;
	else if (option == "output_frequency")
		output_frequency = atoi(value.c_str());
	else if (option == "generations")
		generations = atoi(value.c_str());
	else if (option == "constant_tree")
		constant_tree = atoi(value.c_str());
	else if (option == "tree_type")
		tree_type = atoi(value.c_str());
	else if (option == "max_segment_length")
		max_segment_length = atof(value.c_str());
	else if (option == "substitution_model_type") {
		std::cout <<"Found a model type " << value << std::endl;
		substitution_model_types.push(atoi(value.c_str()));
	} else if (option == "constant_substitution_model") {
		constant_substitution_models.push(atoi(value.c_str()));
	} else if (option == "substitution_model_initialization") {
		substitution_models_initialization_type.push(atoi(value.c_str()));
	} else if (option == "substitution_model_initialization_file") {
		substitution_models_initialization_files.push(value);
	} else if (option == "number_of_substitution_models_in_mixture_model")
		number_of_substitution_models_in_mixture_model = atoi(value.c_str());
	else if (option == "output_directory") {
		output_directory = value;
	} else {
		//STP: Unrecognized option is a non-fatal error so simply print warning
		std::cerr << "Unrecognized option found in control file:" << std::endl;
		std::cerr << "\"" << option << "\" with a value \"" << value << "\""
				<< std::endl;
	}
}

// Might remove this function because it is never called
void Options::WriteOptions(std::ofstream& out_stream) {
	/* Modified 2013-4-1 STP
	 * Using setw(20) cannot handle longer arguments well. In particular, the
	 * output files are long strings and are being cut at 20 characters.
	 * I changed setw to using \t instead
	 */
	out_stream << "Using options specified by control file:" << std::endl;
	out_stream << std::endl << "# debug and loudness controls #" << std::endl;
	out_stream << "\tdebug\t\t\t" << debug << std::endl;
	out_stream << "\tseed\t\t\t" << random_number_generator_seed << std::endl;

	out_stream << std::endl << "# control generations, frequency of sampling #"
			<< std::endl;
	out_stream << "\tgenerations\t\t\t" << generations << std::endl;

	out_stream << std::endl
			<< "# options for simulation to get proposal step in MCMC #"
			<< std::endl;

	out_stream << std::endl << "# probablity calculation method #" << std::endl;

	out_stream << std::endl << "# optional input file names #" << std::endl;
	out_stream << "\ttreefile\t\t\t" << tree_file << std::endl;
	out_stream << "\tseqfile\t\t\t" << sequences_file << std::endl;

	out_stream << std::endl << "# optional output file names #" << std::endl;
	out_stream << "\ttreeoutfile\t\t\t" << tree_out_file << std::endl; // Validated
	out_stream << "\tsuboutfile\t\t\t" << substitutions_out_file << std::endl; // Validated
	out_stream << "\tseqoutfile\t\t\t" << sequences_out_file << std::endl;
	out_stream << "\tlikelihoodoutfile\t\t\t" << likelihood_out_file
			<< std::endl;
}

void Options::UpdateOutputOptions() {
	if (output_directory == "") {
		time_t rawtime = time(NULL); // initializes rawtime
		struct tm * timeinfo = localtime(&rawtime); // transforms it into a local tm struct

		char dir[100];
		strftime(dir, 100, "%Y-%m-%d_%H%M%S/", timeinfo);
		output_directory = std::string(dir);
	}
	char lastchar = output_directory.at(output_directory.length() - 1);
	if (debug)
		std::cout << "The last char in " << output_directory << " is "
				<< lastchar << std::endl;

	if (lastchar != '/' && lastchar != '\\') {
		if (debug)
			std::cout << "last char is not /" << std::endl;
		output_directory += '/';
	}
	// For Windows
	// According to http://sourceforge.net/p/predef/wiki/OperatingSystems/
	// TODO: double check this works on windows for directory that already exists
#ifdef _WIN32
//	if (CreateDirectory(outputDirectory.c_str(), 0)) {
	if (mkdir(output_directory.c_str())) {
		std::cout << "Could not make output directory " << output_directory
				<< std::endl;
	} else {
		std::cout << "Making directory successful" << std::endl;
	}
#else	//ifdef __linux__ || __APPLE__
	struct stat st;
	if (stat(output_directory.c_str(), &st) == 0 && S_ISDIR(st.st_mode)) {
		if (not debug) {
			std::cout << "output dir " << output_directory
			<< " exists. Overwrite? (Y/n)" << std::endl;
			if (getchar() != 'Y')
			exit(1);
		}
	} else {
		mkdir(output_directory.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
	}
	/*Read, write, and search, or execute, for the file owner; S_IRWXU is the bitwise inclusive OR of S_IRUSR, S_IWUSR, and S_IXUSR.
	 Read, write, and search or execute permission for the file's group. S_IRWXG is the bitwise inclusive OR of S_IRGRP, S_IWGRP, and S_IXGRP.
	 Read, write, and search or execute permission for users other than the file owner. S_IRWXO is the bitwise inclusive OR of S_IROTH, S_IWOTH, and S_IXOTH.*/
//		std::cout << "making directory successful: " << m_strOutputDirectory << std::endl;
#endif

	PrependOutputDirectory(tree_out_file);
	PrependOutputDirectory(substitutions_out_file);
	PrependOutputDirectory(likelihood_out_file);
	PrependOutputDirectory(sequences_out_file);
}

void Options::PrependOutputDirectory(std::string &parameter) {
	if (parameter == "")
		std::cerr << "Cannot prepend output directory to empty parameter"
				<< std::endl;
	parameter = output_directory + parameter;
}

void Options::InitializeRandomNumberGeneratorSeed() {
	if (random_number_generator_seed == 0) {
		// If seed not specified, set from clock which is guaranteed to be
		// unique unless multiple runs start at the same time...
		random_number_generator_seed = time(0);
		if (debug) {
			std::cout << "The seed is not specified in the control file"
					<< std::endl;
			std::cout << "Now the seed is " << random_number_generator_seed
					<< std::endl;
		}
	}
	srand(random_number_generator_seed);
}

void Options::CopyFile(string source_filename, string destination_filename) {
	std::ifstream source(source_filename.c_str());
	std::ofstream destination(destination_filename.c_str());

	destination << source.rdbuf();
}

#include "Environment.h"

#include <cstdlib> // For atoi()
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <typeinfo> //Just for testing.

#ifdef _WIN32
#include <dir.h> // For mkdir()
#include <ctime> // For time()
//#include <Windows.h> // For CreateDirectory. Can't use because of CopyFile
#else
#include <sys/stat.h>// For making directories in Linux and OS X
#include <sys/times.h> // For time()
#endif

Environment::Environment() {
	/*
	 *  Initialize environment class.
	 *  Sets a few minor defaults states, however the majority of the default states are 
	 *  read from default.ctrl.
	 */

	seed = 0;
	total_options = 0;
}

// Initialize options.
void Environment::ReadOptions(std::ifstream &default_file_stream, std::ifstream &options_file_stream) {
	ReadControlFile(default_file_stream);
	ReadControlFile(options_file_stream);
	
	debug = get_int("debug");
	InitializeRandomNumberGeneratorSeed();

	PrintOptions();
}

// Read control files.
void Environment::ReadControlFile(std::ifstream &file_stream) {
	std::string key = "";
	std::string value = "";

	std::string line;

	while(std::getline(file_stream, line)) {
		if(line != "") {
			std::istringstream iss(line);
			iss >> key;
			if(key != "#") {
				iss >> value;
				SetOption(key, value);
			}
		}
	}
}

void Environment::SetOption(string option, string value) {
	if(option_to_index.find(option) == option_to_index.end()) {
		//Option does not already exist in options_to_index.
		option_to_index[option] = total_options;
		option_values[total_options] = value;
		total_options++;
	} else {
		//Option does already exist.
		int i = option_to_index[option];
		option_values[i] = value;
	}
}

void Environment::InitializeRandomNumberGeneratorSeed() {
    if (option_to_index.find("seed") != option_to_index.end()) {
	    seed = get_int("seed");
    } else {
	    seed = time(0) ; 
    }
    srand(seed);
}

// Print options.
void Environment::PrintOptions() {

	std::cout << std::endl << "Options:" << std::endl;

	for(std::map<std::string, int>::iterator it = option_to_index.begin(); it != option_to_index.end(); ++it) {
		std::cout << it->first << " = " << option_values[it->second] << " " << std::endl;
	}

	std::cout << std::endl;
}

// Retreiving options.
void Environment::check_option_exists(std::string option) {
	if(option_to_index.find(option) == option_to_index.end()) {
		std::cerr << "Error: Option \"" << option << "\" not set." << std::endl;
		exit(EXIT_FAILURE);
	}
}

std::string Environment::get(std::string option) {
	check_option_exists(option);

	int i = option_to_index[option];
	return option_values[i];
}

int Environment::get_int(std::string option) {
	return atoi(get(option).c_str());
}

float Environment::get_float(std::string option) {
	return atof(get(option).c_str());
}

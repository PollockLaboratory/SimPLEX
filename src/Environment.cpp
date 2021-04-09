#include "Environment.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <time.h>

#include <sys/stat.h> // For making directories in Linux and OS X
#include <sys/times.h> // For time()

Environment::Environment() {
  /*
   *  Initialize environment class.
   *  Sets a few minor defaults states, however the majority of the default states are 
   *  read from default.ctrl.
   */
  seed = 0;
}

// Initialize options.
void Environment::ReadOptions(int argc, char* argv[]) {
  if(argc != 2) {
    std::cerr << "Error: incorrect number of command-line arguments given." << std::endl;
    exit(EXIT_FAILURE);
  }

  std::cout << "Reading options file from: " << argv[1] << std::endl;
  ReadTOMLfile(argv[1]);

  debug = get<bool>("debug");
  ancestral_sequences = get<bool>("INPUT.ancestral_sequences");

  InitializeRandomNumberGeneratorSeed();
}

void Environment::ReadTOMLfile(std::string file_name) {
  try {
    config = cpptoml::parse_file(file_name);
  } catch(const cpptoml::parse_exception& e) {
    std::cerr << "Error: failed to parse TOML configuration file: " << e.what() << std::endl;
    exit(EXIT_FAILURE);
  }
}

void Environment::InitializeRandomNumberGeneratorSeed() {
  seed = get<int>("seed");

  if(seed == 0) {
    seed = time(0) ;
  }

  srand(seed);
}

// Print options.
void Environment::PrintOptions() {
  std::cout << "Received options:" << std::endl << *config << std::endl;
}

// Logger
void Environment::log(std::string message) {
  time_t raw_time;
  time(&raw_time);
  struct tm * timeinfo;
  timeinfo = localtime(&raw_time);
  char buffer[80];
  strftime(buffer, 80, "%D %T", timeinfo);

  //log_stream << "[ " << buffer << " ] " << message << std::endl;
}

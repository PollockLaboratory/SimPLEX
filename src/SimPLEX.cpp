/*
 * SimPLEX
 *
 * Author: Hamish Pike
 * hamish.pike@ucdenver.edu
 * 2023-05-17
 */

#include "SimPLEX.h"
#include "Environment.h"
#include "IO/Files.h"
#include "Model.h"
#include "Data.h"
#include "MCMC.h"

#include "IO/SubstitutionModelParser.h"
#include "ModelParts/SubstitutionModels/SubstitutionModel.h"

#include <sys/times.h>

//Globals
Environment env;
IO::Files files;

double Random() {
  // Returns random number between 0.0 and 1.0;
  return (std::rand() % 10000) / 10000.0;
}

void print_info() {
  std::cout << "SimPLEX - phylogenetic analysis of complex substitution models" << std::endl << std::endl
            << "Usage: SimPLEX configuration_file" << std::endl
            << "\t-v, --version\tshow version of SimPLEX" << std::endl
            << "\t-h, --help\tshow help and overview of SimPLEX" << std::endl;
}

void parse_arguments(int argc, char* argv[]) {
  // an extremely basic argument parser - does nothing except look for version or help.

  if(argc <= 1) {
    print_info();
    std::cout << std::endl << "Error: expecting TOML configuration file." << std::endl;

    exit(EXIT_FAILURE);
  };

  if (strcmp(argv[1], "--version") == 0 or strcmp(argv[1], "-v") == 0) {
    std::cout << "SimPLEX - phylogenetic analysis of complex substitution models" << std::endl
              << "v0.1.0" << std::endl;
    exit(EXIT_SUCCESS);
  } else if (strcmp(argv[1], "--help") == 0 or strcmp(argv[1], "-h") == 0) {
    print_info();

    std::cout << std::endl << "description of SimPLEX" << std::endl;

    std::cout << std::endl << "developed by Hamish NC Pike" << std::endl
              << "hamish.pike@cuanschutz.edu" << std::endl;

    exit(EXIT_SUCCESS);
  } else {
    // assumes the only argument is the configuration file - continue program.
    return;
  }
}

//Entry point for SimPLEX.
int main(int argc, char* argv[]) {
  parse_arguments(argc, argv);
  
  time_t start_time = time(NULL);
  std::cout.precision(17);

  //Establish environment and files.
  env.ReadOptions(argc, argv);

  files.initialize(argv);
  
  std::cout << std::endl;

  // Initiating program.
  Data data;
  data.Initialize();

  std::cout << std::endl;
  
  Model model;
  model.Initialize(data.raw_tree, data.raw_sm);

  data.Uninitialize();

  MCMC mcmc;
  mcmc.Initialize(&model);

  mcmc.Run();

  files.clean_and_close();

  time_t time_taken = time(NULL) - start_time;
  int h = time_taken / 3600;
  int m = (time_taken % 3600) / 60;
  int s = time_taken - (time_taken / 60) * 60;

  char result[100];
  if (h != 0) {
    sprintf(result, "HH:MM:SS %d:%02d:%02d", h, m, s);
  } else {
    sprintf(result, "MM:SS %02d:%02d", m, s);
  }

  std::cout << "Time taken: " << result << std::endl;
  std::string time_taken_file = "Time_taken";

  std::cout << "Successful end." << std::endl;

  return(0);
}


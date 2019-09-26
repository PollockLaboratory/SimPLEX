/*
 * Molecular Evolution Model Testing
 * This is for internal use only.
 *
 * Renamed from SimPLEX on 2014-2-13
 *
 * Author: Hamish Pike
 * hamish.pike@ucdenver.edu
 * 2013-12-05
 */

#include "SimPLEX.h"
#include "Environment.h"
#include "IO/Files.h"
#include "Model.h"
#include "Data.h"
#include "MCMC.h"

#include "IO/SubstitutionModelParser.h"
#include "ModelParts/SubstitutionModels/SubstitutionModel.h"

#include "sol2/sol.hpp"

#include <sys/times.h>

//Globals
Environment env;
IO::Files files;

double Random() {
  // Returns random number between 0.0 and 1.0;
  return (std::rand() % 10000) / 10000.0;
}

//Entry point for SimPLEX.
int main(int argc, char* argv[]) {
  time_t start_time = time(NULL);
  std::cout.precision(17);

  std::cout << std::endl << "SimPLEX" << std::endl
	    << "by Hamish N.C. Pike" << std::endl
	    << "hamish.pike@cuanschutz.edu" << std::endl
	    << "For internal use only." << std::endl
	    << std::endl;

  //Establish environment and files.
  env.ReadOptions(argc, argv);

  files.initialize(argv);
  
  files.add_file("log", env.get<std::string>("OUTPUT.log_out_file"), IOtype::OUTPUT);
  // env.log_stream = files.get_ofstream("log");

  std::cout << "Init complete." << std::endl;

  // Initiating program.
  Data data;
  data.Initialize();

  std::cout << "Data complete." << std::endl;

  Model model;
  model.Initialize(data.raw_tree, data.raw_msa, data.raw_sm);

  std::cout << "Model complete." << std::endl;

  MCMC mcmc;
  mcmc.initialize(&model);

  mcmc.Run();

  model.Terminate();

  files.close();

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


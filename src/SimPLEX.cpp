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
#include "IO.h"
#include "Model.h"
#include "Data.h"
#include "MCMC.h"
#include "utils.h"

#ifdef _WIN32
#include <sys/time.h>
#else
#include <sys/times.h>
#endif

//Globals
Environment env;
IO::Files files;

double Random() {
  return (std::rand() % 10000) / 10000.0;
}

//Entry point for SimPLEX.
int main(int argc, char* argv[]) {
  time_t start_time = time(NULL);
  std::cout.precision(17);

  utils::printHeader();

  //Establish environment and files.
  if(argc == 2) files.set_options_file(argv);
  files.initialize();
  std::ifstream default_file_stream = files.get_ifstream("default");
  std::ifstream options_file_stream = files.get_ifstream("options");
  env.ReadOptions(default_file_stream, options_file_stream);

  files.setupOutputDirectory();

  // Initiating program.
  Data data;
  data.Initialize();

  Model model;
  model.Initialize(data.raw_tree, data.MSA);

  MCMC mcmc;
  mcmc.initialize(&model);

  mcmc.Run();

  model.Terminate();
  utils::Terminate(start_time);

  files.close();

  std::cout << "Successful end." << std::endl;

  return(0);
}


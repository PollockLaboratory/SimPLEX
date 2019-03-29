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

#include "SubstitutionModels/SubstitutionModelParser.h"
#include "SubstitutionModels/SubstitutionModel.h"

#include "sol2/sol.hpp"

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
  env.ReadOptions(argc, argv);

  files.setupOutputDirectory();

  // Initiating program.
  IO::raw_substitution_model* raw_sm = IO::read_substitution_model(env.get<std::string>("DATA.substitution_model_file"));
  SubstitutionModel* sm = new SubstitutionModel();
  sm->from_raw_model(raw_sm);
  sm->finalize();
  
  Data data;
  data.Initialize(sm->get_states());

  Model model;
  model.Initialize(data.raw_tree, data.MSA, sm);

  MCMC mcmc;
  mcmc.initialize(&model);

  mcmc.Run();

  model.Terminate();
  utils::Terminate(start_time);

  files.close();

  std::cout << "Successful end." << std::endl;

  return(0);
}


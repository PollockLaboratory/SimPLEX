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

#include "SubstitutionModels/SubstitutionModelTypes.h"
#include "SubstitutionModels/SubstitutionModel.h"

#include <boost/random/uniform_01.hpp>
#include <boost/random/mersenne_twister.hpp>

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
  //boost::mt19937 rng(time(NULL));
  //static boost::random::uniform_01<boost::mt19937> dist(rng);
  //return(dist());
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
  SubstitutionModel* sm = GetSubstitutionModel();
  sm->Initialize();
  
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


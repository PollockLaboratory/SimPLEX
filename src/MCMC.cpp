#include "MCMC.h"
#include <sstream>

#include "Model.h"
#include "IO/Files.h"

extern double Random();
extern Environment env;
extern IO::Files files;

/// Public Functions ///

MCMC::MCMC() {
  model = 0;
  gen = 0;
  gens = 0;
  lnL = 0;
  complete_likelihood_update = 0;
} 

void MCMC::Initialize(Model* model) {
  /*
   * Init MCMC with model, gens calculate lnL.
   */

  std::cout << "\nInitializing MCMC." << std::endl;
  this->model = model; // associate the pointer with the MCMC

  // Env settings.
  gens = env.get<int>("MCMC.generations");
  out_freq = env.get<int>("MCMC.output_frequency");
  print_freq = env.get<int>("MCMC.print_frequency");
  complete_likelihood_update = env.get<int>("MCMC.full_update_freq");

  //Calculate initial likelihood.
  lnL = model->CalculateLikelihood();

  //Initialize output file.
  files.add_file("likelihoods", env.get<std::string>("OUTPUT.likelihood_out_file"), IOtype::OUTPUT);
  files.write_to_file("likelihoods", "I,GEN,LogL\n");

  RecordState();

  model->printParameters();
  std::cout << "Model and data successfully loaded - MCMC ready." << std::endl << std::endl;
}

void MCMC::sample() {
  sample_status s = model->sample();

  static int gens_since_complete = complete_likelihood_update;

  if(s.full_recalculation or gens_since_complete % complete_likelihood_update == 0) {
    newLnL = model->CalculateLikelihood();
    gens_since_complete = 1;
    //std::cout << "Inter Likelihood: " << newLnL << " ";
  } else {
    // TODO check partial update.

    double delta_LogL = model->CalculateChangeInLikelihood();
    newLnL = lnL + delta_LogL;
    //double test_LogL = model->CalculateLikelihood();
    gens_since_complete++;

    //std::cout << "Inter Likelihood: " << newLnL << " - " << test_LogL << " ";
  }

  if(s.testp) {
    //Metropolis-Hasting method.
    if (log(Random()) <= (newLnL - lnL)) {
      lnL = newLnL;
      model->accept();
      //std::cout << "Accept" << std::endl;
    } else {
      model->reject();
      //std::cout << "Reject" << std::endl;
    }
  } else {
    // No Metropolis Hastings needed - Gibbs sampling.
    lnL = newLnL;
    model->accept();
  }
}

void MCMC::Run() {
  /*
   * Run an initialized MCMC.
   */

  std::cout << "Starting MCMC:" << std::endl;
  for (gen = 1; gen <= gens; gen++) {
    sample();

    if(isnan(lnL)) {
      std::cerr << "Error: LogLikelihood is Nan." << std::endl;
      exit(EXIT_FAILURE);
    }

    if(gen % print_freq == 0) {
      std::cout << "Likelihood: " << lnL << std::endl;
      // model->printParameters();
    }

    if(gen % out_freq == 0) {
      RecordState();
    }
  }
  std::cout << "MCMC complete." << std::endl << std::endl;
}

void MCMC::RecordState() {
  static int i = -1;
  i++;

  std::ostringstream buffer;
  buffer << i << "," << gen << "," << lnL << std::endl;
  files.write_to_file("likelihoods", buffer.str());

  model->RecordState(gen, lnL);
}

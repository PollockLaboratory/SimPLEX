#ifndef MCMC_h_
#define MCMC_h_

#include "Environment.h"
#include "Model.h"

class MCMC {
private:
  Model* model;
  int gen;
  double lnL;
  double newLnL;
  bool accepted;

  // Settings.
  int gens; // Number of generations.
  int out_freq; // Frequency of saving state.
  int print_freq; // Frequency of printing likelihood to term.
  int complete_likelihood_update; // Number of generations before force full likelihood calculation.

  void sample();

  void RecordState(); 
public:
  MCMC();
  void Initialize(Model* model);
  void Run();
};

#endif

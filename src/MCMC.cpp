#include "MCMC.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <limits>
#include "Options.h"
#include "Model.h"

extern double Random();
std::ofstream MCMC::lnlout;


/// Public Functions ///

// default constructor
MCMC::MCMC() {
	model = 0; // should be  model = NULL; ?
	gen = 0;
	gens = 0;
	lnL = 0;
} 
// Init MCMC with model, gens, calculate lnL
void MCMC::Init(Model* model) { std::cout << "Initializing MCMC" << std::endl;
	this->model = model; // associate the pointer with the MCMC
	gens = options.gens;
	lnL = model->CalcLnl(); // should be in Run()?
	lnlout.open(options.lnlout.c_str());
    // Print_Headers(); would be better if this were somewhere else consistent, maybe static vector
	lnlout << "Generation\tLog_likelihood\tProposed_log_likelihood\tAccepted" << std::endl;
}

// Run an initialized MCMC
void MCMC::Run() {  std::cout << "Running MCMC" << std::endl;
	for (gen = 1; gen <= gens; gen++) {
		Model pmodel = ProposeModel(); // returns model, which is prop_mod; bad: newly constructed
		newLnL = pmodel.CalcLnl();

		accepted = TestAccept(newLnL);
		RecordState();

		if (accepted) {  // Danger Will Robinson!
			*model = pmodel; // swap() would be better
			lnL = newLnL;
		}
	}
}

///  Private Functions  ///

// this should be fixed pmodel as parameter of mcmc
Model MCMC::ProposeModel() {
	Model pmodel(*model); // copy current to proposed
	pmodel.SampleParameters(); // propose new parameters
	return pmodel;
}

bool MCMC::TestAccept(double newLnL) {
	return log(Random()) < (newLnL - lnL);
}

void MCMC::RecordState() {  // ought to use a function to return tab separated items with endl
	lnlout << gen << "\t" << lnL << "\t"
		<< newLnL << "\t" << accepted << std::endl;
	model->RecordState(); // is this always getting lnlout?
}

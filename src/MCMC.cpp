#include "MCMC.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <limits>

#include "Model.h"
#include "IO.h"

std::ofstream MCMC::lnlout;

extern double Random();
extern Environment env;
extern IO::Files files;

/// Public Functions ///

MCMC::MCMC() {
	/*
	 * Default constructor.
	 */
	model = 0;
	gen = 0;
	gens = 0;
	lnL = 0;
} 

void MCMC::Init(Model* model) {
	/*
	 * Init MCMC with model, gens calculate lnL.
	 */

	std::cout << "Initializing MCMC" << std::endl;
	this->model = model; // associate the pointer with the MCMC
	std::cout << "Model pointer initialized." << std::endl;
	gens = env.get_int("generations");
	std::cout << "Gens: " << gens << std::endl;

	//Calculate initial likelihood.
	lnL = model->CalcLnl();

	//Initialize output file.
	files.add_file("likelihoods", env.get("likelihood_out_file"), IOtype::OUTPUT);
	lnlout = files.get_ofstream("likelihoods");
	lnlout << "Generation\tLog_likelihood\tProposed_log_likelihood\tAccepted" << std::endl;
}

void MCMC::Run() {
	/*
	 * Run an initialized MCMC.
	 */

	int outfreq = env.get_int("output_frequency");

	std::cout << "Running MCMC" << std::endl;
	for (gen = 1; gen <= gens; gen++) {
		model->SampleParameters();
		newLnL = model->CalcLnl();

		accepted = log(Random()) < (newLnL - lnL);
		if (accepted) { 
			lnL = newLnL;
			model->accept();
		} else {
			model->reject();
		}
		
		if(gen % outfreq == 0) {
			RecordState();
		}
	}
}

///  Private Functions  ///
void MCMC::RecordState() {  // ought to use a function to return tab separated items with endl
	lnlout << gen << "\t" << lnL << "\t" << newLnL << "\t" << accepted << std::endl;
	model->RecordState(); // is this always getting lnlout?
}

#include "MCMC.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <limits>
#include <chrono>

#include "Model.h"
#include "IO.h"

std::ofstream MCMC::lnlout;

extern double Random();
extern Environment env;
extern IO::Files files;

ofstream MCMC::time_out;

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
	tree_sample_freq = env.get_int("tree_sample_frequency");
	
	std::cout << "Gens: " << gens << std::endl;

	//Calculate initial likelihood.
	lnL = model->CalcLnl();
	RecordState();

	//Initialize output file.
	files.add_file("likelihoods", env.get("likelihood_out_file"), IOtype::OUTPUT);
	lnlout = files.get_ofstream("likelihoods");
	lnlout << "I,GEN,LogL" << std::endl;

	files.add_file("time", env.get("time_out_file"), IOtype::OUTPUT);
	time_out = files.get_ofstream("time");

	n_tree_samples = 0;
	total_time_tree_samples = 0;
	n_parameter_samples = 0;
	total_time_parameter_samples = 0;
}

void MCMC::sample() {
	static int i = 1;
	bool sampleType;

	int time_taken;
	std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();

	if(i % tree_sample_freq == 0) {
		float oldLnl = lnL;
		sampleType = model->SampleTree(); // All tree sampling right now is Gibbs.
		lnL = model->CalcLnl();
		//std::cout << "Tree sample: Old: " << oldLnl << " New: " << lnL << " " << (lnL > oldLnl) << std::endl;
		i = 0;

		time_taken = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();
		n_tree_samples++;
		total_time_tree_samples += time_taken;
	} else {
		sampleType = model->SampleSubstitutionModel();
		newLnL = model->CalcLnl();
		if(sampleType) {	
			//Metropolis-Hasting method.
			float r = log(Random());
			//std::cout << "Random: " << r << " " << (newLnL - lnL) << std::endl;
			accepted = r < (newLnL - lnL);
			if (accepted) { 
				lnL = newLnL;
				model->accept();
			} else {
				model->reject();
			}
		} else {
			lnL = newLnL;
		}

		time_taken = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();
		n_parameter_samples++;
		total_time_parameter_samples += time_taken;
	}
	i++;
}

void MCMC::Run() {
	/*
	 * Run an initialized MCMC.
	 */

	int outfreq = env.get_int("output_frequency");
	int print_freq = env.get_int("print_frequency");

	std::cout << "Running MCMC" << std::endl;
	for (gen = 1; gen <= gens; gen++) {
		if(isnan(lnL)) {
			exit(-1);
		}
		if(gen % print_freq == 0) {
			std::cout << "Likelihood: " << lnL << std::endl;
		}
		
		sample();	
		
		if(gen % outfreq == 0) {
			RecordState();
		}
	}

	time_out << "SAMPLING TIME DATA" << std::endl;
	time_out << "Tree Sampling: Total time: " << total_time_tree_samples << " us. Average time per sample: " << total_time_tree_samples/n_tree_samples << " us." << std::endl;
	time_out << "Parameter Sampling: Total time: " << total_time_parameter_samples << " us. Average time per sample: " << total_time_parameter_samples/n_parameter_samples << " us." << std::endl;
	float r = ((float)total_time_tree_samples)/(total_time_tree_samples + total_time_parameter_samples);
	time_out << r*100.0 << "\% of time spent sampling tree parameters." << std::endl;
}

///  Private Functions  ///
void MCMC::RecordState() {  // ought to use a function to return tab separated items with endl
	static int i = -1;
	i++;
	lnlout << i << "," << gen << "," << lnL << std::endl;
	model->RecordState(gen, lnL); // is this always getting lnlout?
}

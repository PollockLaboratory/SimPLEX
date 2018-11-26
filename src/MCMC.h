#ifndef MCMC_h_
#define MCMC_h_

#include <fstream>
#include "Environment.h"
#include "Model.h"

class MCMC {
	public:
		MCMC();
		void Init(Model* model);
		void Run();
	private:
		Model* model;
		int gen;
		double lnL;
		double newLnL;
		bool accepted;

		double newLnL_test; // Temp value for partial likelihood update.

		// Settings.
		int out_freq; // Frequency of saving state.
		int print_freq; // Frequency of printing likelihood to term.
		int gens; // Number of generations.
		int tree_sample_freq; // Frequncy of tree ancestral state sampling.
		
		static std::ofstream lnlout;

		void sample();
		
		void RecordState();
		//bool TestAccept(double newLnL);
		
		// Data for time calculations
		int n_tree_samples;
		int total_time_tree_samples;
		int n_parameter_samples;
		int total_time_parameter_samples;

		static std::ofstream time_out;
};

#endif

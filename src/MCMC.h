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
		
		// Settings.
		int gens;
		int tree_sample_freq;
		
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

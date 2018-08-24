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
		int gens;
		double lnL;
		double newLnL;
		bool accepted;
		
		static std::ofstream lnlout;
		
		void RecordState();
		bool TestAccept(double newLnL);
};

#endif

#ifndef MCMC_h_
#define MCMC_h_

#include <fstream>
#include "Options.h"
#include "Model.h"

extern Options options; // why is this here, exactly

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

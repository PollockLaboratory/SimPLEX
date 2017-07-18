
#ifndef MCMC_h_
#define MCMC_h_

#include <fstream>


#include "Options.h"

//Data is incorporated into the model
//#include "Data.h"
#include "Model.h"

extern Options options;


class MCMC {
public:
	MCMC();

	//Data is incorporated into the model
//	void Initialize(Data* data, Model* model);
	void Initialize(Model* model);
	void Run();

private:
	Model* model;

	//Data is incorporated into the model
	//Data* data;

	double log_likelihood;
	double proposed_log_likelihood;
	bool accepted;

	int generation;
	int generations;

	static std::ofstream likelihood_out;

	void RecordState();

	Model ProposeModel();
	bool IsProposedModelAccepted(double proposed_likelihood);
};

#endif

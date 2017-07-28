
#include "MCMC.h"

#include <iostream>
#include <fstream>
#include <cmath>
#include <limits>

#include "Options.h"
//Data is incorporated into the model
//#include "Data.h"
#include "Model.h"

extern double Random();

std::ofstream MCMC::likelihood_out;

MCMC::MCMC() {
//	data = 0;
	model = 0;
	generation = 0;
	generations = 0;
	log_likelihood = 0;
}

void MCMC::Initialize(Model* model) {
	std::cout << "Initializing MCMC" << std::endl;
	//Data is incorporated into the model
	//this->data = data;
	this->model = model;
	generation = 0;
	generations = options.generations;

	log_likelihood = model->CalculateLogLikelihood();
	
	likelihood_out.open(options.likelihood_out_file.c_str());
	likelihood_out << "Generation\tLog_likelihood\tProposed_log_likelihood\tAccepted" << std::endl;
}

void MCMC::Run() {
	std::cout << "Running MCMC" << std::endl;
	for (generation = 1; generation <= generations; generation++) {
		/* std::cout << "Generation " << generation << " out of " << generations
				<< std::endl; */

		Model proposed_model = ProposeModel();
		proposed_log_likelihood = proposed_model.CalculateLogLikelihood();
		
		accepted = IsProposedModelAccepted(proposed_log_likelihood);
		RecordState();
		
		if (accepted) {
//		if (true) {
			*model = proposed_model;
			log_likelihood = proposed_log_likelihood;
		}
	}
}

Model MCMC::ProposeModel() {
	Model proposed_model(*model);
	proposed_model.SampleParameters();
	return proposed_model;
}

bool MCMC::IsProposedModelAccepted(double proposed_log_likelihood) {
	return log(Random()) < (proposed_log_likelihood - log_likelihood);
}

void MCMC::RecordState() {
	likelihood_out << generation << "\t" << log_likelihood << "\t" 
		<< proposed_log_likelihood << "\t" << accepted << std::endl;
	model->RecordState();
}

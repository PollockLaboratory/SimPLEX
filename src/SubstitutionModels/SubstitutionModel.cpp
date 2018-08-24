#include <iostream>

#include "Environment.h"
#include "IO.h"
#include "SubstitutionModel.h"

extern Environment env;
extern IO::Files files;

SubstitutionModel::SubstitutionModel() {
	substitution_model_out = 0;
}

SubstitutionModel::~SubstitutionModel() {
}

void SubstitutionModel::SampleParameters() {
	/*
	 * Samples a single parameter within the parameter set();
	 */
	parameters.sample();
}

void SubstitutionModel::accept() {
	/*
	 * Accepts the newly sampled parameter set.
	 */
	parameters.accept();
}

void SubstitutionModel::reject() {
	/*
	 * Rejects the newly sampled parameter set, and undoes the changes from the previous sampling.
	 */
	parameters.reject();
}

void SubstitutionModel::RecordState() {
	parameters.RecordStateToFile();
}

void SubstitutionModel::Terminate() {
	delete substitution_model_out;
}

//std::ofstream SubstitutionModel::CreateOutputStream(std::string file_name) {
//	files.add_file("substitution_model", file_name, IOtype::OUTPUT);
//	std::ofstream stream = files.out_file("substitution_model");
//	return(stream);
//}


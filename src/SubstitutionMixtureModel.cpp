#include "SubstitutionMixtureModel.h"

#include <cstdlib> // For exit
#include <iostream>
#include <cmath> // For floor
#include "SubstitutionModelTypes.h"

#include "Options.h"

extern double Random();
extern Options options;

SubstitutionMixtureModel::SubstitutionMixtureModel() {
	// Do I need to initialize any pointers in the substitution_models vector?
	// No
}

SubstitutionMixtureModel* SubstitutionMixtureModel::Clone() {
//	std::cout << "cloning substitution mixture model model" << std::endl;
	return new SubstitutionMixtureModel(*this);
}

SubstitutionMixtureModel::SubstitutionMixtureModel(
		const SubstitutionMixtureModel& substitution_mixture_model) {
//	cout << "SubstitutionMixtureModel copy constructor" << endl;

	substitution_model_out = substitution_mixture_model.substitution_model_out;
	site_assignments = substitution_mixture_model.site_assignments;

	substitution_models.resize(
			substitution_mixture_model.substitution_models.size());
	for (int substitution_model = 0;
			substitution_model < substitution_models.size();
			substitution_model++) {
		substitution_models.at(substitution_model) =
				substitution_mixture_model.substitution_models.at(
						substitution_model)->Clone();
	}
}

//Copy-swap
SubstitutionMixtureModel& SubstitutionMixtureModel::operator=(
		SubstitutionMixtureModel substitution_mixture_model) {
	std::cout << "operator = sub mix mod called " << std::endl;
	std::swap(substitution_model_out,
			substitution_mixture_model.substitution_model_out);
	std::swap(site_assignments, substitution_mixture_model.site_assignments);
	std::swap(substitution_models,
			substitution_mixture_model.substitution_models);

	return (*this);
}

SubstitutionMixtureModel::~SubstitutionMixtureModel() {
	for (int substitution_model = 0;
			substitution_model < substitution_models.size();
			substitution_model++) {
		delete substitution_models.at(substitution_model);
	}
}

void SubstitutionMixtureModel::Initialize(int number_of_sites,
		vector<string> states) {
	SubstitutionModel::Initialize();

	InitializeSiteAssignments(number_of_sites);

	InitializeSubstitutionModels(number_of_sites, states);

	InitializeMixtureModelOut();
}

void SubstitutionMixtureModel::InitializeSiteAssignments(int number_of_sites) {
	site_assignments.resize(number_of_sites);
	SampleSiteAssignments();
}

void SubstitutionMixtureModel::SampleSiteAssignments() {
	for (int site = 0; site < site_assignments.size(); site++) {
		site_assignments.at(site) = std::floor(
				substitution_models.size() * Random());
	}
}

void SubstitutionMixtureModel::InitializeSubstitutionModels(int number_of_sites,
		vector<string> states) {

	substitution_models.resize(
			options.mixture_classes);

	for (int substitution_model = 0;
			substitution_model < substitution_models.size();
			substitution_model++) {
		substitution_models.at(substitution_model) = GetNextSubstitutionModel();
		substitution_models.at(substitution_model)->Initialize(number_of_sites,
				states);
	}
}

void SubstitutionMixtureModel::InitializeMixtureModelOut() {
	string substitution_model_out_file = "substitution_mixture_model";
	options.PrependOutputDirectory(substitution_model_out_file);
	substitution_model_out = new std::ofstream(
			substitution_model_out_file.c_str());

	for (int site = 0; site < site_assignments.size(); site++) {
		if (site != 0)
			*substitution_model_out << "\t";

		*substitution_model_out << "Site_" << site;
	}
	*substitution_model_out << std::endl;
}

void SubstitutionMixtureModel::SampleParameters() {

	SampleSiteAssignments();

	for (int substitution_model = 0;
			substitution_model < substitution_models.size();
			substitution_model++) {
		substitution_models.at(substitution_model)->SampleParameters();
	}
}

void SubstitutionMixtureModel::RecordState() {

	RecordSiteAssignments();

	for (int substitution_model = 0;
			substitution_model < substitution_models.size();
			substitution_model++) {
		substitution_models.at(substitution_model)->RecordState();
	}
}

void SubstitutionMixtureModel::RecordSiteAssignments() {
	for (int site = 0; site < site_assignments.size(); site++) {
		if (site != 0)
			*substitution_model_out << "\t";

		*substitution_model_out << site_assignments.at(site);
	}
	*substitution_model_out << std::endl;
}

double SubstitutionMixtureModel::SubstitutionProbability(int ancestral_state,
		int descendent_state, int site, double branch_length) {

	int site_assignment = site_assignments.at(site);

	double probability =
			substitution_models.at(site_assignment)->SubstitutionProbability(
					ancestral_state, descendent_state, site, branch_length);

	return probability;
}

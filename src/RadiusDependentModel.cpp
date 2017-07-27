#include "RadiusDependentModel.h"

#include <cstdlib> // For exit
#include <iostream>
#include <cmath> // For floor

#include "SubstitutionModelTypes.h"

#include "Options.h"

extern double Random();
extern Options options;

RadiusDependentModel::RadiusDependentModel() {
	// Do I need to initialize any pointers in the substitution_models vector?
	// No
}

RadiusDependentModel* RadiusDependentModel::Clone() {
//	std::cout << "cloning substitution mixture model model" << std::endl;
	return new RadiusDependentModel(*this);
}

RadiusDependentModel::RadiusDependentModel(
		const RadiusDependentModel& substitution_mixture_model) {
//	cout << "RadiusDependentModel copy constructor" << endl;

	substitution_model_out = substitution_mixture_model.substitution_model_out;
	site_radii = substitution_mixture_model.site_radii;

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
RadiusDependentModel& RadiusDependentModel::operator=(
		RadiusDependentModel substitution_mixture_model) {
	std::cout << "operator = sub mix mod called " << std::endl;
	std::swap(substitution_model_out,
			substitution_mixture_model.substitution_model_out);
	std::swap(site_radii, substitution_mixture_model.site_radii);
	std::swap(substitution_models,
			substitution_mixture_model.substitution_models);

	return (*this);
}

RadiusDependentModel::~RadiusDependentModel() {
	for (int substitution_model = 0;
			substitution_model < substitution_models.size();
			substitution_model++) {
		delete substitution_models.at(substitution_model);
	}
}

void RadiusDependentModel::Initialize(int number_of_sites,
		vector<string> states) {
	SubstitutionModel::Initialize();

	site_radii.resize(number_of_sites);
	InitializeState();

	InitializeSubstitutionModels(number_of_sites, states);

	InitializeMixtureModelOut();
}

// These need to be read in the site radii from a file. This should be
// InitializeStateFromFile
void RadiusDependentModel::InitializeStateFromFile(std::string state_in_file) {
	std::ifstream state_in(state_in_file.c_str());

	std::string line;
	getline(state_in, line); // Burn header line

	for (int site = 0; site < site_radii.size(); site++) {
		state_in >> site_radii.at(site);
	}
}

// I might not allow sampling the site radius. This might simply be read in
// from a file and kept constant.
void RadiusDependentModel::SampleSiteRadii() {
	for (int site = 0; site < site_radii.size(); site++) {
		site_radii.at(site) = Random();
	}
}

void RadiusDependentModel::InitializeSubstitutionModels(int number_of_sites,
		vector<string> states) {
	// For now, the number of substitution models in a radius dependent model
	// is fixed at 2: an inner and an outer
	substitution_models.resize(2);

	for (int substitution_model = 0;
			substitution_model < substitution_models.size();
			substitution_model++) {
		substitution_models.at(substitution_model) = GetNextSubstitutionModel();
		substitution_models.at(substitution_model)->Initialize(number_of_sites,
				states);
	}
}

// This will not be needed if the radii are constant. I guess it will be
// needed once at the beginning to record what the initial settings were.
void RadiusDependentModel::InitializeMixtureModelOut() {
	substitution_model_out = new std::ofstream(
			"debug/substitution_mixture_model");

	for (int site = 0; site < site_radii.size(); site++) {
		if (site != 0)
			*substitution_model_out << "\t";

		*substitution_model_out << "Site_" << site;
	}
	*substitution_model_out << std::endl;
}

void RadiusDependentModel::SampleParameters() {

	SampleSiteRadii();

	for (int substitution_model = 0;
			substitution_model < substitution_models.size();
			substitution_model++) {
		substitution_models.at(substitution_model)->SampleParameters();
	}
}

void RadiusDependentModel::RecordState() {

	RecordSiteAssignments();

	for (int substitution_model = 0;
			substitution_model < substitution_models.size();
			substitution_model++) {
		substitution_models.at(substitution_model)->RecordState();
	}
}

// If the site radii are constant then this function might not be required
void RadiusDependentModel::RecordSiteAssignments() {
	for (int site = 0; site < site_radii.size(); site++) {
		if (site != 0)
			*substitution_model_out << "\t";

		*substitution_model_out << site_radii.at(site);
	}
	*substitution_model_out << std::endl;
}

double RadiusDependentModel::SubstitutionProbability(int ancestral_state,
		int descendent_state, int site, double branch_length) {

	double site_radius = site_radii.at(site);

	// substitution_models.at(0) is the inner model and
	// substitution_models.at(1) is the outer model
	double probability = (1 - site_radius)
			* substitution_models.at(0)->SubstitutionProbability(
					ancestral_state, descendent_state, site, branch_length)
			+ site_radius
					* substitution_models.at(1)->SubstitutionProbability(
							ancestral_state, descendent_state, site,
							branch_length);

	return probability;
}

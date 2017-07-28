#ifndef RadiusDependentModel_h_
#define RadiusDependentModel_h_

#include <fstream>
#include <map>
#include <vector>

#include "SubstitutionModel.h"

/**
 * The substitution process depends on the radius as predetermined by the
 * structure and given by the user.
 *
 * This is incomplete.
 *
 * The idea is to begin to incorporate structure information into the
 * evolutionary model. We can determine the radius of a site/amino acid from
 * the structuere and use that information here.
 *
 * How were we going to use it?
 * Perhaps as a prior for whether the site should be assigned to the
 * hydrophobic core model (inner/conserved) or to the hydrophilic model
 * (outer/ fast evolving).
 *
 */

class RadiusDependentModel: public SubstitutionModel {

public:
	RadiusDependentModel();
	RadiusDependentModel(const RadiusDependentModel&);
	RadiusDependentModel& operator=(
			RadiusDependentModel substitution_mixture_model);

	virtual ~RadiusDependentModel();

	virtual RadiusDependentModel* Clone();

	virtual void Initialize(int number_of_sites,
			vector<string> states);

	virtual void SampleParameters();

	virtual void RecordState();

	virtual double SubstitutionProbability(int ancestral_state,
			int descendent_state, int site, double branch_length);
private:
	std::vector<double> site_radii;

	//These must be pointers to take advantage of polymorphism
	// Since these are pointers, I must implement the big three
	std::vector<SubstitutionModel*> substitution_models;

	virtual void InitializeStateFromFile(std::string state_in_file);

	void SampleSiteRadii();
	void InitializeSubstitutionModels(int number_of_sites,
			vector<string> states);
	void RecordSiteAssignments();
	void InitializeMixtureModelOut();
};

#endif

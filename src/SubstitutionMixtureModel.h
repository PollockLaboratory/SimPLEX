#ifndef SubstitutionMixtureModel_h_
#define SubstitutionMixtureModel_h_

#include <fstream>
#include <map>
#include <vector>


#include "SubstitutionModel.h"

/**
 * A mixture model of substitutions models
 *
 */

class SubstitutionMixtureModel: public SubstitutionModel {

public:
	SubstitutionMixtureModel();
	SubstitutionMixtureModel(const SubstitutionMixtureModel&);
	SubstitutionMixtureModel& operator=(SubstitutionMixtureModel substitution_mixture_model);

	virtual ~SubstitutionMixtureModel();

	virtual SubstitutionMixtureModel* Clone();

	virtual void Initialize(int number_of_sites,
					std::vector<std::string> states);

	virtual void SampleParameters();

	virtual void RecordState();

	virtual double SubstitutionProbability(int ancestral_state,
			int descendent_state, int site, double branch_length);
private:
	std::vector<int> site_assignments;

	//These must be pointers to take advantage of polymorphism
	// Since these are pointers, I must implement the big three
	std::vector<SubstitutionModel*> substitution_models;


	void InitializeSiteAssignments(int number_of_sites);
	void SampleSiteAssignments();
	void InitializeSubstitutionModels(int number_of_sites, vector<string> states);
	void RecordSiteAssignments();
	void InitializeMixtureModelOut();

	virtual void InitializeStateFromFile(std::string state_in_file) {};

};

#endif

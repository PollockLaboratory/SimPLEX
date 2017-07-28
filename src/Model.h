/*
 * PLEX specific model for the MCMC. This is separate from the data but knows
 * enough about the data to initialize itself when passed a set of data.
 *
 * This is not a pure model, though, because it will incorporate the data into
 * it. The data in this case is the sequence information at the leaves of the
 * tree. The data is fixed and assumed to be known.
 *
 */

#ifndef Model_h_
#define Model_h_

#include <vector>

#include "Tree.h"

#include "SubstitutionModel.h"

class Model {
public:
	Model();
	Model(const Model& model);
	Model& operator=(Model model);
	~Model();

	void Initialize(map<string, vector<int> > taxa_names_to_sequences,
			vector<string> states);

	void SampleParameters();  // SampleParameters();

	double CalcLnl(); // Calculate the log likelihood

	void RecordState();

	void Terminate();

private:
	static int num_models; // number_of_models;
	int id;

	//Must be a pointer to use polymorphism
	Tree* tree;

	//Must be a pointer to use polymorphism
	SubstitutionModel* sub_model;

	double CalculateLogLikelihoodOfSubtree(Tree& tree);
	double CalculateLogLikelihoodOfChild(Tree& tree, Tree& child);
	void InitSubModel(int number_of_sites,
			vector<string> states);
};

#endif

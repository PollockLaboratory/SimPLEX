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

#include "Trees/Types/Tree.h"
#include "Data.h"

#include "SubstitutionModels/SubstitutionModel.h"

class Model {
public:
	Model();
	~Model();
	void Initialize(IO::RawTreeNode* &raw_tree, SequenceAlignment* &MSA);

	void SampleParameters();
	void accept();
	void reject();
	void printParameters();

	double CalcLnl();
	void RecordState();
	void Terminate();

private:
	Tree* tree;
	SubstitutionModel* substitution_model; //Must be a pointer to use polymorphism

	SubstitutionModel* InitializeSubstitutionModel(int number_of_sites, vector<string> states);
};

#endif

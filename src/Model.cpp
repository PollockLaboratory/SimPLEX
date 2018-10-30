#include <iostream>
#include <cmath>
#include <sys/times.h>
#include <chrono>

#include "Model.h"
#include "Trees/TreeTypes.h"
#include "Trees/TreeParser.h"
#include "SubstitutionModels/SubstitutionModelTypes.h"
#include "SubstitutionModels/SubstitutionModel.h"

#include "Environment.h"
#include "IO.h"

extern Environment env;
extern IO::Files files;

/*  The tree needs the names_to_sequence map. The substitution model might need the empirical frequencies. */

// using namespace std;

Model::Model() {
	/*
	 * The default contructor.
	 */
	tree = NULL;
	SubstitutionModel* substitution_model = NULL;
}

Model::~Model() {
	/*
	 * The destructor function.
	 * Note: we should not be destructing manually.
	 */
	delete tree;
	delete substitution_model;
}

SubstitutionModel* Model::InitializeSubstitutionModel(int num_sites, vector<string> states) {
	/*
	 * This is in the substituionModelType.h, can we make this more explicit?
	 */
	std::cout << "Initializing Substitution Model: ";
	SubstitutionModel* substitution_model = GetSubstitutionModel(); // In SubstituionModelTypes.h
	substitution_model->Initialize(num_sites, states);
	std::cout << std::endl;
	return(substitution_model);
}

void Model::Initialize(IO::RawTreeNode* &raw_tree, SequenceAlignment* &MSA) {
	/*
	 * Initialize the model class.
	 * There are two main components within the model class:
	 * - the tree class - contains the tree topology as well as the sequences.
	 * - the substitution model class - which contains all the rate matrices.
	 */
	int num_sites = (MSA->taxa_names_to_sequences).begin()->second.size();
	substitution_model = InitializeSubstitutionModel(num_sites, MSA->states);
	num_parameters = substitution_model->getNumberOfParameters();

	tree = TreeTypes::pickTreeType();
	tree->Initialize(raw_tree, MSA, substitution_model);
}

bool Model::SampleTree() {
	return(tree->SampleParameters());
}

bool Model::SampleSubstitutionModel() {
	return(substitution_model->SampleParameters());
}

void Model::accept() {
	substitution_model->accept();
}

void Model::reject() {
	substitution_model->reject();
}

void Model::printParameters() {
	tree->printParameters();
}

void Model::RecordState(int gen, double l) {
	/*
	 * Records the state of both the tree and the substitution model.
	 */
	tree->RecordState(gen, l);
	substitution_model->saveToFile(gen, l);
}

double Model::CalcLnl() {
	/*
	 * Calculates the likelihood of the current tree and substitution model.
	 */
	return(tree->calculate_likelihood());
	//return CalculateLogLikelihoodOfSubtree(*tree);
}

void Model::Terminate() {
	/*
	 * Just terminates substitution_model.
	 */
	substitution_model->Terminate();
}

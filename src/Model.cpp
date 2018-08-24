#include "Model.h"

#include <iostream>
#include <cmath>

#include "Trees/TreeTypes.h"
#include "Trees/TreeParser.h"
#include "SubstitutionModels/SubstitutionModelTypes.h"
#include "SubstitutionModels/SubstitutionModel.h"

#include "Environment.h"
#include "IO.h"

extern Environment env;
extern IO::Files files;

/*  The tree needs the names_to_sequence map. The substitution model might need the empirical frequencies. */

using namespace std;

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
	std::cout << "Start initialize." << std::endl;
	SubstitutionModel* substitution_model = GetSubstitutionModel(); // In SubstituionModelTypes.h
	std::cout << "Pointer: " << substitution_model << std::endl;
	substitution_model->Initialize(num_sites, states);
	substitution_model->RecordState();
	std::cout << "Successfully initialized model. " << substitution_model << std::endl;
	return(substitution_model);
}

void Model::Initialize(map<string, vector<int> > taxa_names_to_sequences, vector<string> states) {
	/*
	 * Initialize the model class.
	 * There are two main components within the model class:
	 * - the tree class - contains the tree topology as well as the sequences.
	 * - the substitution model class - which contains all the rate matrices.
	 */
	
	std::string treefile = env.get("tree_file"); 	
	files.add_file("tree_input", treefile, IOtype::INPUT);
	ifstream tree_in = files.get_ifstream("tree_input");

	string tree_string;
	getline(tree_in, tree_string);

	IO::RawTreeNode* raw_tree = IO::parseTree(tree_string);
	tree = TreeTypes::pickTreeType();
	tree->Initialize(raw_tree, taxa_names_to_sequences, states);

	//Not quite sure what this is doing.
	int num_sites = taxa_names_to_sequences.begin()->second.size();
	substitution_model = InitializeSubstitutionModel(num_sites, states);
}

void Model::SampleParameters() {
	/*
	 * Samples both the parameters associated with the tree as well as the substitution model.
	 */
	tree->SampleParameters();
	substitution_model->SampleParameters();
}

void Model::accept() {
	substitution_model->accept();
}

void Model::reject() {
	substitution_model->reject();
}

void Model::RecordState() {
	/*
	 * Records the state of both the tree and the substitution model.
	 */
	tree->RecordState();
	substitution_model->RecordState();
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

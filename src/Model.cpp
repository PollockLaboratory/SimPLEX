#include "Model.h"

#include <iostream>
#include <cmath>

#include "TreeTypes.h"
#include "SubstitutionModelTypes.h"

#include "Options.h"

extern Options options;

/*
 * STP: Exactly what does the model need to be initialized? It will need to
 * initialize the tree and the substitution model. The tree needs the
 * names_to_sequence map. The substitution model might need the empirical
 * frequencies.
 */

using namespace std;

int Model::number_of_models = 0;

Model::Model() {
//	cout << "Model default constructor" << endl;
	id = number_of_models;
	number_of_models++;

	tree = NULL;
	substitution_model = NULL;
}

Model::Model(const Model& model) {
//	cout << "Model copy constructor" << endl;
	id = number_of_models;
	number_of_models++;
	
	// STP: Must include a virtual clone method for derived classes to return 
	// a clone of themselves.
	tree = model.tree->Clone();
	substitution_model = model.substitution_model->Clone();
}

//Copy-swap
Model& Model::operator=(Model model) {
	std::swap(id, model.id);
	std::swap(tree, model.tree);
	std::swap(substitution_model, model.substitution_model);

	return (*this);
}

Model::~Model() {
//	cout << "Model destructor" << endl;
	delete tree;
	delete substitution_model;
}

void Model::Initialize(map<string, vector<int> > taxa_names_to_sequences,
		vector<string> states) {

	tree = InstantiateTree();
	tree->Initialize(taxa_names_to_sequences, states);

	int number_of_sites = taxa_names_to_sequences.begin()->second.size();
	InitializeSubstitutionModel(number_of_sites, states);
}

void Model::InitializeSubstitutionModel(int number_of_sites,
		vector<string> states) {

	substitution_model = GetNextSubstitutionModel();
	substitution_model->Initialize(number_of_sites, states);
}

/**
 *  What if I want the model to be constant? There's not much use for a model
 *  to be constant in an MCMC context. But I might want to use the MCMC to
 *  simply calculate the likelihood. or rather, I might want to use SimPLEX
 *  to just calculate the likelihood of a tree and sequences given a model,
 *  including a substitution model and a tree with the ancestors inferred.
 */

/**
 * Any sort of interaction between the tree and substitution model sampling
 * would happen here. For example, if one wants to use Gibbs sampling for the
 * ancestral states, one would need to know about the substitution model when
 * sampling the states. The only class that knows about both the substitution
 * model and the tree states is the Model class.
 */
void Model::SampleParameters() {
//	std::cout << "Sampling model parameters" << std::endl;
	tree->SampleParameters();
	substitution_model->SampleParameters();
}

void Model::RecordState() {
	tree->RecordState();
	substitution_model->RecordState();
}

double Model::CalculateLogLikelihood() {
	return CalculateLogLikelihoodOfSubtree(*tree);
}

/**
 * This method for calculating the subtree does not result in correct values
 * because the substitution mapping is too simple and incorrect. To properly
 * calculate the likelihood, the substitution model must take into account
 * all the possible substitution mappings that would begin with the ancestral
 * state and end in the descendant state.
 */
double Model::CalculateLogLikelihoodOfSubtree(Tree& tree) {
	double log_likelihood = 0.0;

	if (tree.IsSubtree()) {
		//This is required because children don't know their parent
		log_likelihood += CalculateLogLikelihoodOfChild(tree, *tree.left);
		log_likelihood += CalculateLogLikelihoodOfChild(tree, *tree.right);

		log_likelihood += CalculateLogLikelihoodOfSubtree(*tree.left);
		log_likelihood += CalculateLogLikelihoodOfSubtree(*tree.right);
	}

	return log_likelihood;
}

double Model::CalculateLogLikelihoodOfChild(Tree& tree, Tree& child) {
	double log_likelihood = 0.0;
	for (int site = 0; site < tree.sequence.size(); site++) {
		log_likelihood += log(substitution_model->SubstitutionProbability(
				tree.sequence.at(site), child.sequence.at(site), site,
				child.distance));
	}
	return log_likelihood;
}

void Model::Terminate() {
//tree->terminate
	substitution_model->Terminate();
}

#include "Model.h"

#include <iostream>
#include <cmath>

#include "TreeTypes.h"
#include "SubstitutionModelTypes.h"

#include "Options.h"

extern Options options;

/*  The tree needs the names_to_sequence map. The substitution model might need the empirical frequencies. */

using namespace std;

int Model::num_models = 0; // what is best way to go through models?

// constructor creates tree, sub_model
Model::Model() {   //	cout << "Model default constructor" << endl;
	id = num_models;
	num_models++;
	tree = NULL;
	sub_model = NULL;
}

// copy; we want to avoid this, which creates new model and clones the structures
// STP: Must include a virtual clone method for derived classes to return a clone of themselves.
Model::Model(const Model& model) {   //	cout << "Model copy constructor" << endl;
	id = num_models;  num_models++;
	tree = model.tree->Clone();  sub_model = model.sub_model->Clone();
}

//Copy-swap
Model& Model::operator=(Model model) {
	std::swap(id, model.id);
	std::swap(tree, model.tree);
	std::swap(sub_model, model.sub_model);

	return (*this);
}

// we should not be destructing
Model::~Model() {  //	cout << "Model destructor" << endl;
	delete tree;  delete sub_model;
}

void Model::Initialize(map<string, vector<int> > taxa_names_to_sequences,
		vector<string> states) {

	tree = InstantiateTree();
	tree->Initialize(taxa_names_to_sequences, states);

	int num_sites = taxa_names_to_sequences.begin()->second.size();
	InitSubModel(num_sites, states);
}

void Model::InitSubModel(int num_sites, vector<string> states) {
	sub_model = GetNextSubstitutionModel(); //This is in the substituionModelType.h, can we make this more explicit?
	sub_model->Initialize(num_sites, states);
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
void Model::SampleParameters() { //	std::cout << "Sampling model parameters" << std::endl;
	tree->SampleParameters();
	sub_model->SampleParameters();
}

void Model::RecordState() {
	tree->RecordState();
	sub_model->RecordState();
}

double Model::CalcLnl() {
	return CalculateLogLikelihoodOfSubtree(*tree);
}

/**
 * This method for calculating the subtree does not result in correct values
 * because the substitution mapping is too simple and incorrect. To properly
 * calculate the likelihood, the substitution model must take into account
 * all the possible substitution mappings that would begin with the ancestral
 * state and end in the descendant state.
 */

// make calculations for branches to children, then add calculation for children subtrees
double Model::CalculateLogLikelihoodOfSubtree(Tree& tree) {
	double Lnl = 0.0;
	if (tree.IsSubtree()) {  // required because children don't know their parent
		Lnl += CalculateLogLikelihoodOfChild(tree, *tree.left);
		Lnl += CalculateLogLikelihoodOfChild(tree, *tree.right);

		Lnl += CalculateLogLikelihoodOfSubtree(*tree.left);
		Lnl += CalculateLogLikelihoodOfSubtree(*tree.right);
	}  return Lnl;
}

// foreach site k, get substitution probability from seq state at k in parent to child, over distance (attached to child)
// this needs to be changed for true plex calculation, and avoid log every time
double Model::CalculateLogLikelihoodOfChild(Tree& tree, Tree& child) {  double Lnl = 0.0;
	for (int site = 0; site < tree.sequence.size(); site++) {
		Lnl += log(sub_model->SubstitutionProbability(
				tree.sequence.at(site), child.sequence.at(site), site,
				child.distance));
	}  return Lnl;
}

// model Terminate() just terminates sub_model
void Model::Terminate() {   //tree->terminate
	sub_model->Terminate();
}

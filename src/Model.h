/*
 * PLEX specific model for the MCMC. This is separate from the data but knows
 * enough about the data to initialize itself when passed a set of data.
 */

#ifndef Model_h_
#define Model_h_

#include "SubstitutionCounts.h"
#include "ModelParts/Trees/Tree.h"
#include "ModelParts/SubstitutionModels/SubstitutionModel.h"
#include "ModelParts/ComponentSet.h"
#include "Data.h"

class Model {
public:
  Model();
  void Initialize(IO::RawTreeNode* &raw_tree, IO::raw_substitution_model* &sm);
  sample_status sample();

  void accept();
  void reject();

  double CalculateLikelihood();
  double CalculateChangeInLikelihood();
  double PartialCalculateLikelihood(const double lnL);

  void RecordState(int gen, double l);
  void print();
  void printParameters();
private:
  Valuable* create_uniformization_constant();
  ComponentSet components;

  bool ready; // Checks whether model is ready to be resampled. If not then the changes made from the previous sampling have not been accepted or rejected.
  int num_parameters;

  SubstitutionCounts counts;
  SubstitutionModel* substitution_model;

  Tree* tree;
  CountsParameter* cp;

  std::list<SequenceAlignmentParameter*> msa_parameters;
};

#endif

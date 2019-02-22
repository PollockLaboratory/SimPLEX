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
#include <map>
#include <fstream>

#include "SubstitutionCounts.h"
#include "Trees/Tree.h"
#include "SubstitutionModels/Components/RateVector.h"
#include "Data.h"

#include "SubstitutionModels/SubstitutionModel.h"

class Model {
 public:
  Model();
  ~Model();
  void Initialize(IO::RawTreeNode* &raw_tree, SequenceAlignment* &MSA, SubstitutionModel* &sm);

  bool SampleTree();
  bool SampleSubstitutionModel();

  void accept();
  void reject();

  double CalculateLikelihood();
  double updateLikelihood();
  double PartialCalculateLikelihood(const double lnL);

  void RecordState(int gen, double l);
  void print();
  void printParameters();
  void Terminate();
 private:
  Tree* tree;
  bool ready; // Checks whether model is ready to be resampled. If not then the changes made from the previous sampling have not been accepted or rejected.
  int num_parameters;
  SubstitutionCounts counts;

  SubstitutionModel* substitution_model;
  SubstitutionModel* InitializeSubstitutionModel(int number_of_sites, vector<string> states);

  float u;
  double logL_waiting; // The likelihood of the waiting times.
  double logL_subs; // The likelihood of the waiting times.
  double logL;
  double delta_logL;
};

#endif

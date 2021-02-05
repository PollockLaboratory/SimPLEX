#ifndef TreeParts_h_
#define TreeParts_h_

#include <string>
#include <map>
#include <vector>
#include <iostream>
#include <math.h>

#include "../../IO/TreeParser.h"
#include "../Sequence.h"
#include "../SubstitutionModels/RateVector.h"
#include "../SubstitutionModels/SubstitutionModel.h"
#include "../../SubstitutionCounts.h"

typedef struct Substitution {
  bool occuredp;
  int anc_state;
  int dec_state;
  RateVector* rate_vector; // The rate vector this substitution occured under.
} Substitution;

class TreeNode;

class BranchSegment {
private:
  std::vector<Substitution> substitutions;
  std::vector<RateVector*> rates; //By site.

  inline void update_rate_vectors();
  void set_new_substitutions(); 
public:
  BranchSegment(float distance);
  ~BranchSegment();

  void Initialize(unsigned int n_columns);

  float distance;
  TreeNode* ancestral;
  TreeNode* decendant;

  const std::vector<Substitution>& get_substitutions();
  double get_rate(int pos, int dec_state);

  // Updating.
  void update();
 
  friend std::ostream& operator<< (std::ostream &out, const BranchSegment &b);
};

class TreeNode {
 public:
  static int unique_id;
  std::string name;
  double distance;
  bool sampledp;

  BranchSegment* up;
  BranchSegment* left;
  BranchSegment* right;

  //float** state_probabilities;
  //bool* gaps;

  std::vector<int>* sequence; // Ptr to the sequence.
  SequenceAlignment* MSA; // The MSA that the sequence is in.

  // Hidden States
  std::vector<std::vector<int>*> hidden_state_sequences;

  SubstitutionModel* SM;

  TreeNode();
  TreeNode(IO::RawTreeNode* raw_tree);
  TreeNode(std::string n);

  void connect_substitution_model(SubstitutionModel*);

  //TreeNode* set_gaps();
  
  //Sampling
  bool ready_to_sample();

   // Printing/Display
  std::string toString();
  std::string get_sequence();
  std::string state_at_pos(int i);

  bool isTip();
};

#endif

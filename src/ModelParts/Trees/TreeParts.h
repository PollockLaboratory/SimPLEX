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

  float** state_probabilities;
  bool* gaps;

  std::vector<int>* sequence;
  SequenceAlignment* MSA;

  SubstitutionModel* SM;

  TreeNode();
  TreeNode(IO::RawTreeNode* raw_tree);
  TreeNode(std::string n);

  TreeNode* set_gaps();
  
  //Sampling
  bool ready_to_sample();
  int pick_state_from_probabilities(int pos);
  void calculate_state_probabilities_pos(int pos, TreeNode* left, TreeNode* right, TreeNode* up);

  TreeNode* calculate_state_probabilities(const std::list<int>&);
  void pick_sequences(const std::list<int>&);

  // Printing/Display
  std::string toString();
  std::string get_sequence();
  std::string state_at_pos(int i);

  bool isTip();
};

#endif

#ifndef TreeParts_h_
#define TreeParts_h_

#include <string>
#include <map>
#include <vector>
#include <iostream>
#include <math.h>

#include "../../IO/TreeParser.h"

//class SequenceAlignment;
class SubstitutionModel;
class RateVector;

typedef struct Substitution {
  bool occuredp;
  int anc_state;
  int dec_state;
  RateVector* rate_vector; // The rate vector this substitution occured under.
} Substitution;

class TreeNode;

class BranchSegment {
private:
  unsigned int n_pos;

  std::map<std::string, std::vector<RateVector*>> rates; // Rate vectors by site also.
  std::map<std::string, std::vector<Substitution>> substitutions; // Temp - should be private.

  inline void update_rate_vectors();
  void set_new_substitutions(); 
public:
  float distance;
  TreeNode* ancestral;
  TreeNode* decendant;

  BranchSegment(float distance);
  ~BranchSegment();

  void Initialize(unsigned int n_columns, std::map<std::string, std::list<std::string>> all_states);
  const std::vector<Substitution>& get_substitutions(std::string domain);

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

  // Hidden States
  std::map<std::string, std::vector<signed char>*> sequences; // Name -> sequence.

  SubstitutionModel* SM;

  TreeNode();
  ~TreeNode();
  TreeNode(IO::RawTreeNode* raw_tree);
  TreeNode(std::string n);

  void connect_substitution_model(SubstitutionModel*);

  // Printing/Display
  std::string toString();
  //std::string get_sequence();
  //std::string state_at_pos(int i);

  unsigned long get_hash_state(unsigned int);
  unsigned long get_hypothetical_hash_state(unsigned int pos, std::string domain, signed char state);

  bool isTip();
};

#endif

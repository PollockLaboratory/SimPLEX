#ifndef TreeParts_h_
#define TreeParts_h_

#include <string>
#include <map>
#include <vector>
#include <iostream>
#include <math.h>

#include "../../IO/TreeParser.h"
#include "../SubstitutionModels/States.h"

//class SequenceAlignment;
class SubstitutionModel;
class RateVector;

typedef struct Substitution {
  bool occuredp;
  state_element anc_state;
  state_element dec_state;
  RateVector* rate_vector; // The rate vector this substitution occured under.
} Substitution;

class TreeNode;

class BranchSegment {
private:
  unsigned int n_pos;

  // String Key is the state domain.
  // TODO should be in ordered containers like a vector rather than a map.
  std::map<std::string, std::vector<RateVector*>> rates; // Rate vectors by site also.
  std::map<std::string, std::vector<Substitution>> substitutions; // Substitutions by site

  inline void update_rate_vectors();
  void set_new_substitutions();

  // Finding applicable RateVectors.
  unsigned long get_hypothetical_hash_state(std::map<std::string, state_element>& states, unsigned int pos);
public:
  class iterator;

  float distance;
  TreeNode* ancestral;
  TreeNode* decendant;

  BranchSegment(float distance);
  ~BranchSegment();

  void Initialize(unsigned int n_columns, std::list<std::string> state_domain_names);

  const Substitution& get_substitution(std::string domain, unsigned int pos);
  const std::vector<Substitution>& get_substitutions(std::string domain);

  // Finding applicable RateVectors.
  RateVector* get_hypothetical_rate_vector(std::string focal_domain, std::map<std::string, state_element>& states, unsigned int pos);

  // Loop through substitutions.
  BranchSegment::iterator begin(unsigned int pos);
  BranchSegment::iterator end();

  // Update susbtitution counts and rate vectors.
  void update();
 
  friend std::ostream& operator<< (std::ostream &out, const BranchSegment &b);

  class iterator:public std::iterator<std::output_iterator_tag, std::pair<std::string, Substitution>> {
  public:
    explicit iterator(BranchSegment&, unsigned int pos, bool end);
    std::pair<std::string, Substitution> operator*() const;
    iterator& operator++(int);
    bool operator!=(const iterator &) const; 
  private:
    BranchSegment& branch; 
    unsigned int pos;
    std::map<std::string, std::vector<Substitution>>:: iterator it;
    std::pair<std::string, Substitution> value;
  };
};

class TreeNode {
 public:
  static int unique_id;
  std::string name;
  double distance;

  // This is a general utility to mark a node temporarily.
  // For example, useful when recursively navigating the tree and marking previously visited nodes.
  bool tagp;

  BranchSegment* up;
  BranchSegment* left;
  BranchSegment* right;

  // State Domain Sequences.
  std::map<std::string, std::vector<state_element>*> sequences; // domain_name -> sequence.

  SubstitutionModel* SM;

  TreeNode();
  ~TreeNode();
  TreeNode(IO::RawTreeNode* raw_tree);
  TreeNode(std::string n);

  void connect_substitution_model(SubstitutionModel*);

  unsigned long get_hash_state(unsigned int);

  // Printing/Display
  std::string toString();


  bool isTip();
};

#endif

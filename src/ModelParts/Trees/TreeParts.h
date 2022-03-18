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
  signed char anc_state;
  signed char dec_state;
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
public:
  class iterator;

  float distance;
  TreeNode* ancestral;
  TreeNode* decendant;

  BranchSegment(float distance);
  ~BranchSegment();

  void Initialize(unsigned int n_columns, std::map<std::string, std::list<std::string>> all_states);

  const Substitution& get_substitution(std::string domain, unsigned int pos);
  const std::vector<Substitution>& get_substitutions(std::string domain);

  // NEW
  signed char get_alt_domain_state(std::string alt_domain, std::string view_domain, unsigned int pos);
  unsigned long get_hypothetical_hash_state(std::string domain, signed char state, unsigned int pos);
  unsigned long get_hypothetical_hash_state(std::string focal_domain, std::map<std::string, signed char>& states, unsigned int pos);

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

  // State Domains
  std::map<std::string, std::vector<signed char>*> sequences; // domain_name -> sequence.

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
  //unsigned long get_hypothetical_hash_state(unsigned int pos, std::string domain, signed char state);

  bool isTip();
};

#endif

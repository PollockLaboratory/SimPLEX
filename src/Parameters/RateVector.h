#ifndef RateVector_h_
#define RateVector_h_

#include <vector>
#include <set>
#include <string>
#include <list>

#include "AbstractValue.h"
#include "ContinuousFloat.h"
#include "VirtualSubstitutionRate.h"

class BranchSegment; // Defined in Trees/Types/TreeParts.h

struct bpos {
	// Represents the position and branch that a rate vector applies to.
	BranchSegment* branch;	
	int pos;

	bool operator<(const bpos& x) const {
		if(branch == x.branch) {
			if(pos < x.pos) {
				return(true);
			} else {
				return(false);
			}
		}
		if(branch < x.branch) return(true);
	}
};


// Fundamental collection type of parameters in substitution model.
class RateVector {
	public: 
		RateVector(std::string, int state, std::vector<AbstractValue*>);

		std::string name;
		std::vector<AbstractValue*> rates;
		int state; // Determines the (ancestral) state that this rate vector applies to.
		
		void remove_location(int pos, BranchSegment* bs);
		void add_location(int pos, BranchSegment* bs);
		std::set<bpos> get_locations();

		void print();
	private:
		std::set<bpos> locations; // List of branches that rate vector applies to.
		int size;
		static int IDc;
};

// Collections of rate vectors.
class RateVectorSet {
	// This collection holds all of the rate vectors currently in the substitution model.
	public:
		RateVectorSet();
		std::vector<RateVector*> c;
		void Initialize();
		
		RateVector*& operator[] (const int i);

		void add(RateVector* v);

		void print();
		void saveToFile(int gen, double l);
	private:
		static std::ofstream out_file;
};

#endif

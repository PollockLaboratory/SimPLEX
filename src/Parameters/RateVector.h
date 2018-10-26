#ifndef RateVector_h_
#define RateVector_h_

#include <vector>
#include <string>

#include "AbstractValue.h"
#include "ContinuousFloat.h"
#include "VirtualSubstitutionRate.h"

// Fundamental collection type of parameters in substitution model.
class RateVector {
	public: 
		RateVector(std::string, int state, std::vector<AbstractValue*>);
		RateVector(std::string name, int size, int state, float u);
		std::vector<AbstractValue*> rates;
		int state; // Determines the state that this rate vector applies to.
		std::string name;
		void print();
	private:
		int size;
		static int IDc;
		inline void create_parameters(int n, float u);
		
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

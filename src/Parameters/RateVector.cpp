#include "RateVector.h"
#include "Environment.h"

extern Environment env; 

int RateVector::IDc = 0;

inline void RateVector::create_parameters(int n, float u) {
	VirtualSubstitutionRate* unifp = new VirtualSubstitutionRate(this->name + "-virtual", u);
	for(int i = 0; i < n; i++) {
		if(i == state) {
			rates.push_back(unifp);
		} else {
			std::string name = this->name + "-" + std::to_string(i);
			ContinuousFloat* p = new ContinuousFloat(name, 0.1, 0.3);
			unifp->add_dependancy(p);
			rates.push_back(p);
		}
	}
	unifp->refresh();
}

// RateVector::RateVector(int size) {
//	this->size = size;
//	name = std::to_string(IDc);
//	IDc++;	
//	create_parameters(size);
//}

RateVector::RateVector(std::string name, int size, int state, float u) {
	/*
	 * The state variable also describes the position of the virtual subtitution rate in the matrix.
	 */
	this->size = size;
	this->state = state;
	this->name = name;
	create_parameters(size, u);
}

void RateVector::print() {
	std::cout << "RateVector:\t" << name << "\t";
	for(auto it = rates.begin(); it != rates.end(); ++it) {
		std::cout << (*it)->getValue() << " ";
	}
	std::cout << std::endl;
}

// COLLECTIONS of rate vectors.
// RateMatrix.

int RateMatrix::IDc = 0;

inline void RateMatrix::create_vectors(int n) {
	float u = env.get_float("uniformization_constant");
	std::cout << "Unif: " << u << std::endl;
	for(int i = 0; i < n; i++) {
		std::cout << "Creating rate vectors." << std::endl;
		std::string name = this->name + "-" + std::to_string(i);
		RateVector* v = new RateVector(name, n, i, u);
		rv.push_back(v);
	}
}

RateMatrix::RateMatrix(int size) {
	this->size = size;
	name = std::to_string(IDc);
	IDc++;
	create_vectors(size);
}

RateMatrix::RateMatrix(std::string name, int size) {
	this->size = size;
	this->name = name;
	create_vectors(size);
}

// RateVectorSet

RateVectorSet::RateVectorSet() {
}

void Initialize() {
}

RateVector*& RateVectorSet::operator[] (const int i) {
	return(c[i]);
}

void RateVectorSet::add(RateVector* v) {
	c.push_back(v);
}

void RateVectorSet::add(RateMatrix* Q) {
	for(std::vector<RateVector*>::iterator it = Q->rv.begin(); it != Q->rv.end(); ++it) {
		c.push_back(*it);
	}
}

void RateVectorSet::print() {
	for(std::vector<RateVector*>::iterator it = c.begin(); it != c.end(); ++it) {
		(*it)->print();
	}
}


#include "RateVector.h"
#include "Environment.h"

extern Environment env; 

int RateVector::IDc = 0;

RateVector::RateVector(std::string name, int state, std::vector<AbstractValue*> params) : name(name), state(state) {
	size = params.size();
	rates = params;
}

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

void RateVectorSet::print() {
	for(std::vector<RateVector*>::iterator it = c.begin(); it != c.end(); ++it) {
		(*it)->print();
	}
}


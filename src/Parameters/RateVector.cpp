#include "RateVector.h"
#include "Environment.h"
#include "IO.h"

extern Environment env;
extern IO::Files files;

int RateVector::IDc = 0;

// Figure out locations Data Structure.
RateVector::RateVector(std::string name, int state, std::vector<AbstractValue*> params) : name(name), state(state) {
	size = params.size();
	rates = params;	
	locations = {};
}

void RateVector::add_location(int pos, BranchSegment* bs) {
	bpos loc = {bs, pos};
	locations.insert(loc);
}

void RateVector::remove_location(int pos, BranchSegment* bs) {
	bpos loc = {bs, pos};
	locations.erase(loc);
}

std::set<bpos> RateVector::get_locations() {
	    return(locations);
}

// Util.
void RateVector::print() {
	std::cout << "RateVector:\t" << name << "\t";
	for(auto it = rates.begin(); it != rates.end(); ++it) {
		std::cout << (*it)->getValue() << " ";
	}
	std::cout << std::endl;
}

// COLLECTIONS of rate vectors.
std::ofstream RateVectorSet::out_file;

RateVectorSet::RateVectorSet() {
}

void RateVectorSet::Initialize() {
	files.add_file("rate_vectors", env.get("rate_vectors_out_file"), IOtype::OUTPUT);
	out_file = files.get_ofstream("rate_vectors");
		
	out_file << "I,GEN,LogL,NAME,ANC";
	for(auto it = env.state_to_integer.begin(); it != env.state_to_integer.end(); ++it) {
		out_file << "," << it->first;
	}
	out_file << std::endl;
}

RateVector*& RateVectorSet::operator[] (const int i) {
	return(c[i]);
}

void RateVectorSet::add(RateVector* v) {
	for(int i = 0; i < v->rates.size(); i++) {
		v->rates[i]->state = i;
		v->rates[i]->rv = v;
	}
	c.push_back(v);
}

void RateVectorSet::print() {
	for(std::vector<RateVector*>::iterator it = c.begin(); it != c.end(); ++it) {
		(*it)->print();
	}
}

void RateVectorSet::saveToFile(int gen, double l) {
	static int i = -1;
	++i;
	for(auto it = c.begin(); it != c.end(); ++it) {
		out_file << i << "," << gen << "," << l << "," << (*it)->name << "," << (*it)->state;
		for(auto jt = (*it)->rates.begin(); jt != (*it)->rates.end(); ++jt) {
			out_file << "," << (*jt)->getValue();
		}
	out_file << std::endl;
	}
}

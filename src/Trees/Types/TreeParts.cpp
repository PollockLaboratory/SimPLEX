/*
 * TREE PARTS
 * This file holds the classes for the components of the tree:
 * 	- Tree nodes - Contains the sequences.
 * 	- Branche segments - Contains the branch lengths and location of substitutions.
 */

#include "TreeParts.h"
#include "Environment.h"
#include <iostream>

extern double Random();
extern Environment env;

// Branch segment.
BranchSegment::BranchSegment(float distance) {
	this->distance = distance;
	rates = {};
//	std::cout << "Making new branch segment. Distance: " << this->distance << std::endl;
}

BranchSegment::~BranchSegment() {
	decendant->up = 0;
}

std::ostream& operator<< (std::ostream &out, const BranchSegment &b) {
	out << b.distance;
	return out;
}

bool BranchSegment::virtualSubstituionQ(int state) {
	float u = env.u;
	double rate = decendant->SM->selectRateVector(state)->rates[state]->getValue();

	double noSub = 1.0/(1.0 + u*distance);
	double Sub = (rate * distance)/(1 + u*distance);

	noSub = noSub/(Sub+noSub);
	Sub = Sub/(Sub+noSub);

	if(Random() < Sub) {
		return(true);
	} else {
		return(false);
	}
}

void BranchSegment::updateStats() {
	num0subs = 0;
	num1subs = 0;
	subs = {};

	std::vector<int> anc = *(ancestral->sequence);
	std::vector<int> dec = *(decendant->sequence);
	substitution s;
	for(int i = 0; i < anc.size(); i++) {
		if(anc.at(i) == dec.at(i)) {
			if(virtualSubstituionQ(anc.at(i))) {
				num1subs += 1;
				s = {i, anc.at(i), dec.at(i)};
				subs.push_back(s);
			} else {
				num0subs += 1;
			}
		} else {
			num1subs += 1;
			s = {i, anc.at(i), dec.at(i)};
			subs.push_back(s);
		}
	}
}

// Tree nodes.

int TreeNode::unique_id = 0;

TreeNode::TreeNode(IO::RawTreeNode* raw_tree) {
	/* 
	 * TreeNode Constructor - from rawTreeNode.
	 */
	name = raw_tree->name;
	distance = raw_tree->distance;
	up = 0;
	left = 0;
	right = 0;
	sampled = false;
}

TreeNode::TreeNode(std::string n) {
	/* 
	 * TreeNode Constructor.
	 */
	name = n;
	up = 0;
	left = 0;
	right = 0;
	sampled = false;
}

TreeNode::TreeNode() {
	/* 
	 * TreeNode Constructor - from name and distance.
	 */
	name = "Node" + std::to_string(unique_id);
	unique_id++;
	up = 0;
	left = 0;
	right = 0;
	sampled = false;
}

// Sampling.
void branchLikelihood(double &l, int anc, int dec, float t_b, SubstitutionModel* SM) {
	/*
	 * Calculates the likelihood of a branch, and adds it to the vector l.
	 */
	// Doesn't take into account virtual substitutions.
	float u = env.u;
	if(anc == dec) {
		l *= 1.0/(1.0 + t_b *u);
	} else {
		double rate = SM->selectRateVector(anc)->rates[dec]->getValue();
		l *= (rate*t_b)/(1.0 + t_b *u);
	}
}

std::vector<double> normalizeLikelihoods(std::vector<double> &l) {
	double total = 0.0;
	for(auto it = l.begin(); it != l.end(); ++it) {
		total += *it;
	}

	for(auto it = l.begin(); it != l.end(); ++it) {
		*it = *it / total;
	}

	return(l);
}

void TreeNode::sampleSinglePosition(int pos) {
	// Don't hard code 20 states
	std::vector<double> l(20, 1.0);	
	for(int state = 0; state < 20; state++) {
		if(up) {
			int anc = up->ancestral->sequence->at(pos);
			branchLikelihood(l[state], anc, state, up->distance, SM);
		}

		if(left) {
			int dec = left->decendant->sequence->at(pos);
			branchLikelihood(l[state], state, dec, left->distance, SM);
		}

		if(right) {
			int dec = right->decendant->sequence->at(pos);
			branchLikelihood(l[state], state, dec, right->distance, SM);
		}
	}
	
	// std::cout << "Un-normalized likelihoods: ";
	//for(auto it = l.begin(); it != l.end(); ++it) {
	 //	std::cout << *it << " ";
	//}
	//std::cout << std::endl;

	//std::cout << "Normalized likelihoods: ";
	l = normalizeLikelihoods(l);
	//for(auto it = l.begin(); it != l.end(); ++it) {
	//	std::cout << *it << " ";
	//}
	
	double r = Random();
	int i = 0;
	double c = l[0];
	while(r > c) {
		i++;
		c += l[i];
	}
	
	//if(sequence->at(pos) != i) {
//		std::cout << "Old aa: " << sequence->at(pos) << " New: " << i << std::endl;
//	}

	(*sequence)[pos] = i;
}

void TreeNode::sampleSequence() {
	if(!isTip()) {
		for(int i = 0; i < sequence->size(); i++) {
			sampleSinglePosition(i);		
		}
	}
}

bool TreeNode::isTip() {
	if(this->left == 0 and this->right == 0) {
		return(true);
	} else {
		return(false);
	}
}



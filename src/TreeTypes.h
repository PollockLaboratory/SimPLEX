#ifndef TreeTypes_h_
#define TreeTypes_h_

#include "Tree.h"
#include "Tree_B1.h"

#include "Options.h"
extern Options options;

/**
 * I'm not sure of the best way to do this. These are not globals; they only
 * are needed in one file: Model.cpp
 *
 */

static const int normal_tree_type = 0;
static const int b1_tree_type = 1;

static Tree* InstantiateTree() {
	Tree* tree = NULL;

	// This is a good candidate for a switch or case construction
	if (options.tree_type == normal_tree_type) {
		tree = new Tree();
	} else if (options.tree_type == b1_tree_type) {
		tree = new Tree_B1();
	} else {
		std::cerr
				<< "Tree type not recognized "
				<< options.tree_type << std::endl;
		std::exit(-1);
	}
	return tree;
}

#endif

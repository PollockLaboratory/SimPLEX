#include <iostream>

#include "TreeTypes.h"

#include "Environment.h"
extern Environment env;

static const int default_tree_type = 0;

Tree* TreeTypes::pickTreeType() {
	/*
	 * Looks in the environment variables to identify the correct tree class.
	 */

	Tree* tree = NULL;
	int chosen_tree_type = env.get_int("tree_type");
	
	std::cout << "INSTANTIATING THE TREE:" << std::endl;
	if (chosen_tree_type == default_tree_type) {
		tree = new Tree();
	} else {
		std::cerr << "Tree type not recognized " << chosen_tree_type << std::endl;
		std::exit(EXIT_FAILURE);
	}
	return tree;
}



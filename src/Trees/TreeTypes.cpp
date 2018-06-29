#include <iostream>

#include "Types/Tree.h"
#include "Types/Tree_B1.h"

#include "Options.h"
extern Options options;

static const int normal_tree_type = 0;
static const int b1_tree_type = 1;

Tree* InstantiateTree() {
	Tree* tree = NULL;
	int chosen_tree_type = options.get_int("tree_type");
	
	std::cout << "INSTANTIATING THE TREE:" << std::endl;
	if (chosen_tree_type == normal_tree_type) {
		tree = new Tree();
	} else if (chosen_tree_type == b1_tree_type) {
		tree = new Tree_B1();
	} else {
		std::cerr
				<< "Tree type not recognized "
				<< chosen_tree_type << std::endl;
		std::exit(-1);
	}
	return tree;
}



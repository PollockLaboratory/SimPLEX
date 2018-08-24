#ifndef TreeTypes_h_
#define TreeTypes_h_

#include "Types/Tree.h"

/**
 * I'm not sure of the best way to do this. These are not globals; they only
 * are needed in one file: Model.cpp
 *
 */
namespace TreeTypes {
	Tree* pickTreeType();
};

#endif

#ifndef UTILS_h_
#define UTILS_h_

#ifdef _WIN32
#include <sys/time.h>
#else
#include <sys/times.h>
#endif

#include <string>

namespace utils {
	void printHeader();
	void Terminate(time_t start_time);
}

#endif

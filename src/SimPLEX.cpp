/*
 * Molecular Evolution Model Testing
 * This is for internal use only.
 *
 * Renamed from SimPLEX on 2014-2-13
 *
 * Author: Stephen Pollard
 * Stephen.T.Pollard@gmail.com
 * 2013-12-05
 */

#include "SimPLEX.h"
#include "Options.h"
#include "Model.h"
#include "Data.h"
#include "MCMC.h"
#include "utils.h"

#ifdef _WIN32
#include <sys/time.h>
#else
#include <sys/times.h>
#endif

//Globals
Options options;

double Random() {
	return (std::rand() % 10000) / 10000.0;
}

//Entry point for SimPLEX.
int main() {
	time_t start_time = time(NULL);

	utils::printHeader();
	options.ReadOptions();

	Data data;
	data.Initialize();

	Model model;
	model.Initialize(data.taxa_names_to_sequences, data.states);

	MCMC mcmc;
	mcmc.Init(&model);

	mcmc.Run();

	// Now that I have objects allocated on the heap, this is required. Unless I use shared pointers...
	model.Terminate();
	utils::Terminate(start_time);

	std::cout << "Successful end" << std::endl;

	return 0;
}


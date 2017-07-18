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



Options options;

int main() {
	SimPLEX::Initialize();

	Data data;
	data.Initialize();

	Model model;
	model.Initialize(data.taxa_names_to_sequences, data.states);

	MCMC mcmc;
	mcmc.Initialize(&model);

	mcmc.Run();

	// Now that I have objects allocated on the heap, this is required. Unless
	// I use shared pointers...
	model.Terminate();

	SimPLEX::Terminate();

	return 0;
}


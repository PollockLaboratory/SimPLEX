#ifndef SubstitutionModel_h_
#define SubstitutionModel_h_

#include <fstream>
#include <sstream>

#include <vector>
using std::vector;
#include <map>
using std::map;
#include <string>
using std::string;

#include "../Parameters/ParameterSet.h"
#include "../Parameters/RateVector.h"

class SubstitutionModel {
	public:
		SubstitutionModel();
		virtual void Initialize(int number_of_sites, std::vector<std::string> states) = 0;

		RateVector* selectRateVector(int state);

		bool SampleParameters();
		void accept(); //After a model is sampled it must be accepted or rejected before next sampling.
		void reject();

		void printParameters();
		int getNumberOfParameters();

		void saveToFile(int gen, double l);
		virtual void Terminate();
	protected:
		std::ofstream* substitution_model_out;

		void add_rate_vector(RateVector* v);

		void finalize();

		//std::ofstream CreateOutputStream(std::string file_name);
	private:
		ParameterSet parameters;
		RateVectorSet rateVectors;
		
};

#endif

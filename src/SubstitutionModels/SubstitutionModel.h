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
#include "TreeParts.h"

class SubstitutionModel {
	public:
		SubstitutionModel();
		virtual void Initialize(int number_of_sites, std::vector<std::string> states) = 0;

		RateVector* selectRateVector(BranchSegment* b, int pos);

		void SampleParameters();
		void accept(); //After a model is sampled it must be accepted or rejected before next sampling.
		void reject();
		void printParameters();

		void RecordState();
		virtual void Terminate();
	protected:
		std::ofstream* substitution_model_out;

		void add_rate_vector(RateVector* v);
		void add_rate_matrix(RateMatrix* Q);

		void finalize();

		//std::ofstream CreateOutputStream(std::string file_name);
	private:
		ParameterSet parameters;
		RateVectorSet rateVectors;
		

};

#endif

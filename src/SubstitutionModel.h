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


/**
 *
 * The purpose of a substitution model is to calculate the probability of
 * a substitution or not a substitution. It must be able to give the
 * probability of not seeing a substitution at a site given a certain branch
 * length and starting amino acid.
 *
 * Perhaps it should be generic and simply take all information needed to
 * calculate the probability.
 *
 *
 * Alternate names: EvolutionModel
 *
 */

class SubstitutionModel {

public:
	SubstitutionModel();
	virtual ~SubstitutionModel();

	// Must provide this if using polymorphism
	virtual SubstitutionModel* Clone() = 0;

	/**
	 * I'm not sure if I will have a single initialize function that takes in
	 * all information that is relevent for the most complex model and then
	 * have the simpler models ignore that information or simply have multiple
	 * Initialize functions with different parameters.
	 *
	 * Either way, I need a way to pass more information to more complex models
	 * that need that information for initialization.
	 *
	 * This is getting out of hand, I think I want a single initialize function
	 * that takes all the information that any of the substitution models
	 * might need. The derived classes can choose to ignore some of the
	 * information.
	 *
	 */

	virtual void Initialize(int number_of_sites,
			std::vector<std::string> states) = 0;

	/**
	 * I have made these functions purely virtual so that derived classes must
	 * implement them. But for a constant class, these functions will be empty.
	 * Do I want constant classes to have to implement these functions even if
	 * they are trivial?
	 *
	 * The trade off is that if they are purely virtual, one must explicitly
	 * defined these function as doing nothing. The explicitness is good but
	 * it adds work. For now I will choose the option that requires more work.
	 * If it ends up being way too much work then I will make them not purely
	 * virtual.
	 */

	virtual void SampleParameters() = 0;

	virtual void RecordState() = 0;

	/**
	 * Similarly to the initialize functions, I'm not sure if I will provide
	 * a bunch of substitution probability functions here all taking different
	 * arguments or simply a single function and certain models ignore certain
	 * arguments.
	 */
	virtual double SubstitutionProbability(int ancestral_state,
			int descendent_state, int site, double branch_length) = 0;

	virtual void Terminate();

protected:
	// Can this be private instead?
	// No because derived classes need access to them

	int id;

	bool is_constant; // Could this functionality be accomplished using
	// constant objects instead and making a SampleParameters() constant
	// method that guarantees the object will not change?

	// I could use a shared pointer here. What would be the advantages?
	std::ofstream* substitution_model_out;

	void Initialize(); // Notice this is not virtual
	/**
	 * InitializeState() is not virtual because all substitution models will use
	 * this function and alter InitializeStateFromFile().
	 */
	void InitializeState(); // Notice this is not virtual
	/**
	 * This should be purely virtual because I want to require derived classes
	 * to implement this. They have not all implemented it yet and so it is not
	 * virtual yet.
	 */
	virtual void InitializeStateFromFile(std::string state_in_file) = 0;

	std::string IdToString();
};

#endif

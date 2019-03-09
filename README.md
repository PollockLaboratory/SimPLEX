# SimPLEX
Simplex is a tool for phylogenetic analysis, focussing on fitting complex substitution models to large datasets.

It takes the input of a rooted tree, sequences for the node tips and a substitution model, and returns the posterior distribution of the parameters inside that substitution model.

## Installation

	$ cd build/

	$ cmake ..

## Usage

When configuring simPLEX there are 3 main core input files.

1. Tree file - a rooted tree file in the newark tree format.
2. Sequence file - a fasta file with the sequences for each of taxa on the tips of the tree.
3. Model file - a lua script that configures the substitution model.
4. Configuration file - a configuration for the Markov Chain parameters.

## Configuration Parameters

* **debug** - bool - when true(1) will print extra error messages. This is not really that widely used at the moment.
* **max_segment_length**- float - the maximum length of a branch segment. The branch splitting algorithm will split branches til all segments are below this length.
* **branch_split_algorithm** - int - options: 0 for no branch splitting (not reccommended), 1 for split in half til under max_segment_length.
* **uniformization_constant** - float - the uniformatization constant, must be between 0 and 1.
* **ancestral_sequences** - bool - if the ancestral sequences have already been calculated. TRUE(1) if they have, OFF(0) if they have not.
* **substitution_model_type** - int - select the substitution model type.
* **custom_model** - file path - location of the lua file of the custom model, will only be read if custom_model is selected through substitution_model_type.
* **tree_sample_frequency** - int - the frequency for the ancestral sequence to be resampled.
* **generations** - int - number of generations for the Markov chain.
* **output_frequency** - int - the frequency at which the state of the Markov chain will be saved to the output files.
* **print_frequency** - int - the frequency at which the log likelihood will be printed to the command line. This is primarily a debug tool/sanity check: you can watch the log likelihood increasing over your chain.
* **tree_file** - file path - the input tree file. Newick tree format.
* **sequence_file** - file path - the input sequences file. Fasta format.
* **tree_out_file** - file name - the name of the file the labelled tree should be outputted to.
* **sequences_out_file** - file name - the name of the file for the ancestrally reconstructed sequences.
* **parameter_out_file** - file name - the name of the file for the rate parameters.
* **likelihoods_out_file** - file name - the name of the file for the likelihoods file.

# Lua model script.

This is where most of the configuration is done. The model is constructed as a set of rate vectors that can apply to certain locations on the tree.

## Functions
* **model**
  * set_name
  * add_rate_vector
* **states**
  * set
  * print
  * to_int
  * to_str
  
## Parameters
* **Continuous Float**
* **Category Float**
* **Fixed Float**
* **Rate Vector**

# Structure of the program

Here is a description of the implimentation of the program. The main intent of this section is to assist anyone who might be interested in editing the code, it is not neccesary for a end user to know this infomation.

## Tree
The tree structure has two type of nodemplimented as two separate classes:

* Tree Nodes - These hold the sequences. Each has 3 pointers: UP going to the branch to the ancestor; LEFT and RIGHT going to the two decendent branchSegments. There are 4 types of TreeNode specified by whether the pointers are connected to anything.
  * Tip - only the UP pointer is not NULL
  * Internal Node - All 3 nodes point to branch segments.
  * Link Internal Node - Only LEFT and UP are connected.
  * Root - Only the LEFT and RIGHT pointers are connected.
* BranchSegements - These have the lengths, substitutions and pointers to the ratevector that apply. They also have pointers to the ancestral and decendant nodes.

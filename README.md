# SimPLEX

Simplex is a tool for phylogenetic analysis, focussing on fitting complex substitution models to large datasets.

It takes the input of a rooted tree, sequences for the node tips and a substitution model, and returns the posterior distribution of the parameters inside that substitution model.

## Installation

### Download a copy of SimPLEX source code

	$ git clone https://github.com/PollockLaboratory/SimPLEX.git

### Dependencies
simPLEX has a number of external dependencies that can sometimes to be challenging to link/compile into the final binary. These include lua5.2 and the C++ boost library. To install them on Ubuntu/Debian based linux distributions use:

	$ sudo apt update
	
	$ sudo apt install liblua5.2-dev libboost-all-dev
	
### Compilation
The simplest way to install this script is to use the install script. This will automatically carry out the following commands and copy the simPLEX binary to the /usr/bin directory on the user PATH.

	$ sudo ./install.sh
	
To install simPLEX manually, use the following commands: 

	$ mkdir build
	
	$ cd build/

Set compiler options.

	$ cmake .. -DCMAKE_BUILD_TYPE=DEBUG	

or, for full optimizations:
	
	$ cmake .. -DCMAKE_BUILD_TYPE=RELEASE

Actually compile:

	$ make
	
Test that everything is working:

	$ ./bin/simPLEX ../resources/options.toml
	
Copy the simPLEX binary to a location on your path so simPLEX can be used in any directory, for example:

	$ cp ./bin/simPLEX /usr/bin/simPLEX
	
## Usage

When configuring simPLEX there are 3 main core input files.

1. Tree file - a rooted tree file in the newark tree format.
2. Sequence file - a fasta file with the sequences for each of taxa on the tips of the tree.
3. Model file - a lua script that configures the substitution model.
4. Configuration file - a configuration for the Markov Chain parameters.

### Configuration Parameters

* **debug** - bool - when true(1) will print extra error messages. This is not really that widely used at the moment.
* **max_segment_length**- float - the maximum length of a branch segment. The branch splitting algorithm will split branches til all segments are below this length.
* **branch_split_algorithm** - int - options: 0 for no branch splitting (not reccommended), 1 for split in half til under max_segment_length.
* **threshold** - the smallest size that a virtual substitution rate can be - related to the uniformization constant.
* **step_size** - the maximum size step the uniformization constant can be.
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

## Citations

Bioinformatics. 2012 Nov 15;28(22):2989-90. doi: 10.1093/bioinformatics/bts555. Epub 2012 Sep 12. Phylogenetics, likelihood, evolution and complexity A P Jason de Koning 1 , Wanjun Gu, Todd A Castoe, David D Pollock, PMID: 22976081 PMCID: PMC3496332 DOI: 10.1093/bioinformatics/bts555

Mol Biol Evol. 2010 Feb;27(2):249-65. Epub 2009 Sep 25. Rapid likelihood analysis on large phylogenies using partial sampling of substitution histories A P Jason de Koning 1 , Wanjun Gu, David D Pollock PMID: 19783593 PMCID: PMC2877550 DOI: 10.1093/molbev/msp228

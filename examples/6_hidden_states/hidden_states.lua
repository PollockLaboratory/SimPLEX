-- Markov state
nucleotide_states = { "a", "t", "c", "g" }
States.new("nucleotide", nucleotide_states,
           {sequences_output = Config.get_str("MODEL.sequences_out"), substitutions_output = Config.get_str("MODEL.substitutions_out")})
Data.load_state("nucleotide", Config.get_str("MODEL.sequences_in"))

hidden_states = { "A", "B" }
States.new("hidden", hidden_states,
	   { sequences_output = Config.get_str("MODEL.hidden_sequence_out"), substitutions_output = Config.get_str("MODEL.hidden_substitutions_out")})
Data.load_state("hidden", Config.get_str("MODEL.hidden_sequences_in"));

-- Rate Vectors for Nucleotide state dimension.
Q_hidden_A = {}
Q_hidden_B = {}

for i=1,#nucleotide_states do
   Q_hidden_A[i] = {}
   Q_hidden_B[i] = {}
end

A_rate = Parameter.new("base-A", "continuous", {initial_value = 0.01, step_size = Config.get_float("MODEL.step_size"), lower_bound = 0.0});
B_rate = Parameter.new("base-B", "continuous", {initial_value = 0.01, step_size = Config.get_float("MODEL.step_size"), lower_bound = 0.0});

for i=1,#nucleotide_states do
	nucleotide_i = States["nucleotide"][i]
	for j=i,#nucleotide_states do
	   if i ~= j then
	      Q_hidden_A[i][j] = A_rate
	      Q_hidden_A[j][i] = A_rate

	      Q_hidden_B[i][j] = B_rate
	      Q_hidden_B[j][i] = B_rate
	   else
	      Q_hidden_A[i][j] = Parameter.new("A-virtual-"..tostring(nucleotide_i), "virtual", {})
	      Q_hidden_B[i][j] = Parameter.new("B-virtual-"..tostring(nucleotide_i), "virtual", {})
	   end
	end

	Model.add_rate_vector(RateVector.new("RV-nucleotide-A-"..nucleotide_i,
					     {domain = "nucleotide", state = nucleotide_i, hidden="A", pos = {}},
					     Q_hidden_A[i]))
	Model.add_rate_vector(RateVector.new("RV-nucleotide-B-"..nucleotide_i,
					     {domain = "nucleotide", state = nucleotide_i, hidden="B", pos = {}},
					     Q_hidden_B[i]))
end

-- Rate Vectors for Hidden state dimension.
hidden_rate = Parameter.new("hidden-rate", "fixed", {value = 0.01})

Model.add_rate_vector(RateVector.new("RV-hidden-A",
				     {domain = "hidden", state = "A", nucleotide = "*", pos = {}},
				     {Parameter.new("virtual-A", "virtual", {}), hidden_rate}))

Model.add_rate_vector(RateVector.new("RV-hidden-B",
				     {domain = "hidden", state = "B", nucleotide = "*", pos = {}},
				     {hidden_rate, Parameter.new("virtual-B", "virtual", {})}))



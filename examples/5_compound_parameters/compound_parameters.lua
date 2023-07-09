--[[

   GTR - General time reversable model
   The substitution rate from amino acid i to j is the product of the equlibrium frequency and the mutability.

   rate_i_j = equil_i * m_i_j
   where:
     rate_i_j = substitution rate from nucleotide i to j.
     equil_i = equilibrium frequency of nucleotide i.
     m_i_j = mutability of nucleotide i to j
]]

Model.set_name("GTR")

model_states = Config.get_string_array("MODEL.states")

States.new("nucleotide", model_states,
	   {sequences_output = Config.get_str("MODEL.sequences_out_file"), substitutions_output = Config.get_str("MODEL.substitutions_out_file")})

Data.load_state("nucleotide", Config.get_str("MODEL.sequences_file"))

-- Create parameters for equilibrium frequencies.
equil_freq = {}
Q = {}
for i=1,#model_states do
   equil_freq[i] = Parameter.new("freq-"..tostring(States.nucleotide[i]), "continuous", {initial_value = 0.001, step_size = Config.get_float("MODEL.step_size"), lower_bound = 0.0})
   Q[i] = {}
end

for i=1,#model_states do
	for j=i,#model_states do
	   if i ~= j then

	      print(i, j)

	      rate = Parameter.new("base-"..tostring(States.nucleotide[i])..tostring(States.nucleotide[j]), "continuous", {initial_value = 0.001, step_size = Config.get_float("MODEL.step_size"), lower_bound = 0.0})
	      Q[i][j] = Parameter.named_multiply(tostring(States.nucleotide[i])..tostring(States.nucleotide[j]), rate, equil_freq[j])
	      Q[j][i] = Parameter.named_multiply(tostring(States.nucleotide[j])..tostring(States.nucleotide[i]), rate, equil_freq[i])
	   else
	      Q[i][j] = Parameter.new("virtual-"..tostring(States.nucleotide[i]), "virtual", {})
	   end
	end
	Model.add_rate_vector(RateVector.new("RV-"..tostring(States.nucleotide[i]), {domain = "nucleotide", state = States.nucleotide[i], pos = {}}, Q[i]))
end

--[[

   True GTR

]]

Model.set_name("GTR")

model_states = Config.get_string_array("MODEL.states")

print(1)

States.new("nucleotide", model_states,
	   {sequences_output = Config.get_str("MODEL.sequences_out_file"), substitutions_output = Config.get_str("MODEL.substitutions_out_file")})

Data.load_state("nucleotide", Config.get_str("MODEL.sequences_file"))

-- Create parameters for equilibrium frequencies.
print(2)

equil_freq = {}
Q = {}
for i=1,#model_states do
   equil_freq[i] = Parameter.new("freq-"..tostring(States.nucleotide[i]), "continuous", {initial_value = 0.001, step_size = Config.get_float("MODEL.step_size"), lower_bound = 0.0})
   Q[i] = {}
end

print(3)

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

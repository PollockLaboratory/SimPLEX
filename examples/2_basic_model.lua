--[[

   (ANOTHER) VERY SIMPLE SUBSTITUTION MODEL - for nucleotide sequences.

   This script builds on the model presented in the previous example, but rather than using a predefined
   list of states, expands the model to any list of n states defined in the configuration file.

   Thus the new general model is:
   q_ij = x, where i != j
   q_ij = -(n-1)x = v_i, where i == j, where v_i is the virtual substitution rate.

]]

Model.set_name("JC69 - general")

model_states = Config.get_string_array("MODEL.states")

States.new("nucleotide", model_states,
	   {sequences_output = Config.get_str("MODEL.sequences_out_file"), substitutions_output = Config.get_str("MODEL.substitutions_out_file")})

Data.load_state("nucleotide", Config.get_str("MODEL.sequences_file"))

x = Parameter.new("x", "continuous", {initial_value = 0.02, step_size = Config.get_float("MODEL.step_size"), lower_bound = 0.0})

for i=1,#model_states do
	rv_list = {}
	for j=1,#model_states do
	   if i ~= j then
	      rv_list[j] = x
	   else
	      rv_list[j] = Parameter.new("virtual-"..tostring(States.nucleotide[i]), "virtual", {})
	   end
	end
	Model.add_rate_vector(RateVector.new("RV-"..tostring(States.nucleotide[i]), {domain = "nucleotide", state = States.nucleotide[i], pos = {}}, rv_list))
end

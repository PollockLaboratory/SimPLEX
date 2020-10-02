--[[

   (ANOTHER) VERY SIMPLE SUBSTITUTION MODEL - for nucleotide sequences.

   This script builds on the model presented in the previous example, but rather than using a predefined
   list of states, expands the model to any list of n states defined in the configuration file.

   Thus the new general model is:
   q_ij = x, where i != j
   q_ij = -(n-1)x = v_i, where i == j, where v_i is the virtual substitution rate.

]]

model.set_name("JC69 - general")

model_states = config.get_string_array("MODEL.states")
states.set(model_states)

x = Parameter.new("x", "continuous", {initial_value = 0.001, step_size = config.get_float("MODEL.step_size"), lower_bound = 0.0})

for i=1,#model_states do
	rv_list = {}
	for j=1,#model_states do
	   if i ~= j then
	      rv_list[j] = x
	   else
	      rv_list[j] = Parameter.new("virtual-"..tostring(model_states[i]), "virtual", {})
	   end
	end
	model.add_rate_vector(RateVector.new("RV-"..tostring(model_states[i]), {state = model_states[i], pos = {}}, rv_list))
end

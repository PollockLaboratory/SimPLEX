model.set_name("GTR")

state_names = config.get_string_array("MODEL.states")

states.set(state_names)

-- Setup Codons.

codon_table = {A = {},
	       R = {},
	       N = {},
	       D = {},
	       C = {},
	       E = {},
	       Q = {},
	       G = {},
	       H = {},
	       I = {},
	       L = {},
	       K = {},
	       M = {},
	       F = {},
	       P = {},
	       S = {},
	       T = {},
	       W = {},
	       Y = {},
	       V = {}}

for k, v in pairs(state_names) do
   print(k, v)
end

param_template = {initial_value = 0.001, step_size = 0.0005, lower_bound = 0}

print(#state_names)

Q = {}

for i=1,#state_names do
	Q[i] = {}
	for j=1,#state_names do
	   Q[i][j] = Parameter.new(tostring(state_names[i])..tostring(state_names[j]), "continuous", param_template)
	end
	model.add_rate_vector(RateVector.new("RV-"..tostring(state_names[i]), {state = states[i], pos = {}}, Q[i]))
end

model.set_name("GTR")

state_names = config.get_string_array("MODEL.states")

states.set(state_names)

for k, v in pairs(state_names) do
   print(k, v)
end

-- Setup Matrix
equil_freqs = {}
Q = {}
for i=1,#state_names do
   equil_freqs[i] = Parameter.new("freq_" .. tostring(states[i]), "continuous", {initial_value = 0.1, step_size = 0.001, lower_bound = 0.0})
   Q[i] = {}
end

print("setup equilibrium frequencies")

for i=1,#state_names do
   for j=i,#state_names do
      if i == j then
	 Q[i][i] = Parameter.new(tostring(states[i])..tostring(states[j]), "fixed", {value = 0.1})
      else
	 mutability = Parameter.new("mut"..tostring(states[i])..tostring(states[j]), "continuous", {initial_value = 0.01, step_size = 0.001, lower_bound = 0.0})
	 Q[i][j] = Parameter.named_multiply(tostring(states[i])..tostring(states[j]), equil_freqs[i], mutability)
	 Q[j][i] = Parameter.named_multiply(tostring(states[j])..tostring(states[i]), equil_freqs[j], mutability)
      end
   end
end

print("Setup rate matrix")

for i=1,#state_names do
   model.add_rate_vector(RateVector.new("RV-"..tostring(state_names[i]), {state = states[i], pos = {}}, Q[i]))
end


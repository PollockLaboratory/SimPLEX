model.set_name("Joseph's Model")

state_names = config.get_string_array("MODEL.states")
states.set(state_names)

for k, v in pairs(state_names) do
   print(k, v)
end

param_template = {initial_value = 0.02}

rate = Parameter.new("rate", "float", {initial_value = 0.005})

Q = {}
for i=1,20 do
	Q[i] = {}
	for j=1,20 do
		-- Q[i][j] = Parameter.new(tostring(i)..tostring(j), param_template)
		if i == j then
		   Q[i][j] = Parameter.new(tostring(state_names[i])..tostring(state_names[j]), "float", param_template)
		else
		   --Q[i][j] = Parameter.new(tostring(state_names[i])..tostring(state_names[j]), "float", {initial_value = 0.005})
		   Q[i][j] = rate
		end
	end
	model.add_rate_vector(RateVector.new("RV-"..tostring(state_names[i]), {state = states[i], pos = {}}, Q[i]))
end

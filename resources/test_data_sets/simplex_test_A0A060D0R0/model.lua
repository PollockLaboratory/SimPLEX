model.set_name("Joseph's Model")

states.set(config.get_string_array("MODEL.states"))

param_template = {initial_value = 0.01}

rate = Parameter.new("rate", "float", {initial_value = 0.005})

Q = {}
for i=1,20 do
	Q[i] = {}
	for j=1,20 do
		-- Q[i][j] = Parameter.new(tostring(i)..tostring(j), param_template)
		if i == j then
			Q[i][j] = Parameter.new(tostring(i)..tostring(j), "float", param_template)
		else
			Q[i][j] = rate
		end
	end
	model.add_rate_vector(RateVector.new("RV-"..tostring(i), {state = states[i], pos = {}}, Q[i]))
end

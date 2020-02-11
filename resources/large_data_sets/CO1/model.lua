model.set_name("Single rate model.")

state_names = config.get_string_array("MODEL.states")
parameter_type = config.get_str("MODEL.parameter_type")

states.set(state_names)

for k, v in pairs(state_names) do
   print(k, v)
end

if parameter_type == "continuous" then
   print("Continuous")
   rate = Parameter.new("rate", "continuous", {initial_value = 0.01, step_size = 0.0001})
elseif parameter_type == "discrete" then
   print("Discrete")
   cats = {}
   for i = 0.0015, 0.003, 0.00015 do
      cats[#cats+1] = i
   end
   
   rate = Parameter.new("rate", "discrete", {categories = cats})
else
   print("Error: no matching parameter type.")
end

param_template = {initial_value = 0.001, step_size = 0.0001}

print(#state_names)

Q = {}
for i=1,#state_names do
	Q[i] = {}
	for j=1,#state_names do
		-- Q[i][j] = Parameter.new(tostring(i)..tostring(j), param_template)
		if i == j then
		   Q[i][j] = Parameter.new(tostring(state_names[i])..tostring(state_names[j]), "continuous", param_template)
		else
		   --Q[i][j] = Parameter.new(tostring(state_names[i])..tostring(state_names[j]), "float", {initial_value = 0.005})
		   Q[i][j] = rate
		end
	end
	model.add_rate_vector(RateVector.new("RV-"..tostring(state_names[i]), {state = states[i], pos = {}}, Q[i]))
end

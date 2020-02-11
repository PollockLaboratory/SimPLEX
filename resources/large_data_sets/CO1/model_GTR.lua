model.set_name("Single rate model.")

state_names = config.get_string_array("MODEL.states")
parameter_type = config.get_str("MODEL.parameter_type")

states.set(state_names)

for k, v in pairs(state_names) do
   print(k, v)
end

if parameter_type == "continuous" then
   print("Continuous")
elseif parameter_type == "discrete" then
   print("Discrete")
   cats = {}
   for i = 0.0015, 0.003, 0.00015 do
      cats[#cats+1] = i
   end   
else
   print("Error: no matching parameter type.")
end

param_template = {initial_value = 0.01, step_size = 0.001}

print(#state_names)

-- Setup Matrix
Q = {}
for i=1,#state_names do
   Q[i] = {}
end

for i=1,#state_names do
   for j=i,#state_names do
      Q[i][j] = Parameter.new(tostring(states[i])..tostring(states[j]), parameter_type, param_template)
      Q[j][i] = Q[i][j]
   end
end

for i=1,#state_names do
   model.add_rate_vector(RateVector.new("RV-"..tostring(state_names[i]), {state = states[i], pos = {}}, Q[i]))
end


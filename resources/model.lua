model.set_name("General Time Reversible")

states.set(config.get_string_array("MODEL.states"))

param_template = {type = "float", initial_value = 0.01}

Q1 = {}
Q2 = {}

for i = 1, states.count do
   Q1[i] = {}
   Q2[i] = {}
   for j = 1, states.count do
      --print(i, j, states[i], states[j])
      if i <= j then
	 Q1[i][j] = Parameter.new(states[i]..states[j]..1, param_template)
	 Q2[i][j] = Parameter.new(states[i]..states[j]..2, param_template)
      else
	 Q1[i][j] = Q1[j][i]
	 Q2[i][j] = Q1[j][i]
      end
   end

   model.add_rate_vector(RateVector.new("RV1-"..states[i], {state = states[i]}, Q1[i]))
   --model.add_rate_vector(RateVector.new("RV2-"..states[i], {state = states[i]}, Q2[i]))

end

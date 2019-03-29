model.set_name("General Time Reversible")

states.set(config.get_string_array("MODEL.states"))

param_template = {type = "float", initial_value = 0.01}

Q = {}

for i = 1, states.count do
   Q[i] = {}
   for j = 1, states.count do
      --print(i, j, states[i], states[j])
      if i <= j then
	 Q[i][j] = Parameter.new(states[i]..states[j], param_template)
      else
	 Q[i][j] = Q[j][i]
      end
   end

   model.add_rate_vector(RateVector.new("RV-"..states[i], {state = states[i]}, Q[i]))

end

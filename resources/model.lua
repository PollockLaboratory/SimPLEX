model.set_name("General Time Reversible")

states.set({"A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"})
states.print()

Q = {}

for i = 1, 20 do
   Q[i] = {}
   for j = 1, 20 do
      if i == j then
	 Q[i][j] = VirtualSubstitutionRate.new(states.to_str(i) .. states.to_str(j))
      elseif i < j then
	 Q[i][j] = ContinuousFloat.new(states.to_str(i) .. states.to_str(j), 0.01, 0.001, 0.0)
	 -- Q[i][j] = CategoryFloat.new(states.to_str(i) .. states.to_str(j), {0.01, 0.02, 0.03, 0.04})
	 -- Q[i][j] = FixedFloat.new(states.to_str(i) .. states.to_str(j), 0.01)
      else
	 Q[i][j] = Q[j][i]
      end
   end

   for j = 1, 20 do
      if i ~= j then
	 Q[i][i]:add_rate(Q[i][j])
      end
   end

   model.add_rate_vector(RateVector.new("RV-"..states.to_str(i), i, Q[i]))

end

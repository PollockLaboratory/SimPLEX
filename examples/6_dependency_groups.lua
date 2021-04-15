--[[

   True GTR

]]

-- This has not been fully implimented at all.

Model.set_name("GTR")

model_states = Config.get_string_array("MODEL.states")
States.set(model_states)

-- Create parameters for equilibrium frequencies.

dgroup = DependencyGroup.new("test", {})

equil_freq = {}
Q = {}
for i=1,#model_states do
   equil_freq[i] = Parameter.new("freq-"..tostring(States[i]), "continuous", {initial_value = 0.001, step_size = Config.get_float("MODEL.step_size"), lower_bound = 0.0})
   Q[i] = {}
end

for i=1,#model_states do
	for j=i,#model_states do
	   if i ~= j then
	      rate = Parameter.new("base-"..tostring(States[i])..tostring(States[j]), "continuous", {initial_value = 0.001, step_size = Config.get_float("MODEL.step_size"), lower_bound = 0.0})
	      Q[i][j] = Parameter.named_multiply(tostring(States[i])..tostring(States[j]), rate, equil_freq[j])
	      Q[j][i] = Parameter.named_multiply(tostring(States[j])..tostring(States[i]), rate, equil_freq[i])
	   else
	      Q[i][j] = Parameter.new("virtual-"..tostring(States[i]), "virtual", {})
	   end
	end
	Model.add_rate_vector(RateVector.new("RV-"..tostring(States[i]), {state = States[i], pos = {}}, Q[i]))
end

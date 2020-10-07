model.set_name("GNR-name")

state_names = config.get_string_array("MODEL.states")

states.set(state_names)

-- Setup Categories

lower_bound = 0.0
upper_bound = 0.003
step_size = 0.0004

cats = {}

for i = 0, (upper_bound - lower_bound)/step_size do
   l = lower_bound + (i * step_size)
   w = lower_bound + ((i + 1) * step_size)
   init = (l + w)/2.0

   cats[i+1] = Parameter.new("Category-" .. i, "continuous", {initial_value = init, step_size = 0.01, lower_bound = 0})
end

Category_Object = Categories.new("RateCategories", cats)
   
Q = {}

for i=1,#state_names do
	Q[i] = {}
	for j=1,#state_names do
	   Q[i][j] = Parameter.new(tostring(state_names[i])..tostring(state_names[j]), "discrete", {categories = Category_Object})
	end
	model.add_rate_vector(RateVector.new("RV-"..tostring(state_names[i]), {state = states[i], pos = {}}, Q[i]))
end

model.set_name("GNR-name")

state_names = config.get_string_array("MODEL.states")

states.set(state_names)

n = config.get_int("MODEL.num_categories")

-- Setup Categories
lower_bound = 0.0
upper_bound = 0.3
step_size = (upper_bound - lower_bound) / n 

cats = {}

for i = 0, n - 1  do
   l = lower_bound + (i * step_size)
   w = lower_bound + ((i + 1) * step_size)
   init = (l + w)/2.0

   cats[i+1] = Parameter.new("Category-" .. i, "continuous", {initial_value = init, step_size = 0.01, lower_bound = 0})
end

for i = 1, #cats do
   if i ~= 1 then
      cats[i]:set_lower_bound(cats[i-1])
   end

   if i ~= #cats then
      cats[i]:set_upper_bound(cats[i+1])
   end
end

Category_Object = Categories.new("RateCategories", cats)

Q = {}

for i=1,#state_names do
	Q[i] = {}
	for j=1,#state_names do
	   Q[i][j] = Parameter.new(tostring(state_names[i])..tostring(state_names[j]), "discrete", {categories = Category_Object})
	   --Q[i][j] = Parameter.new(tostring(state_names[i])..tostring(state_names[j]), "continuous", {initial_value = 0.01, step_size = 0.1, lower_bound = 0.0})
	end
	model.add_rate_vector(RateVector.new("RV-"..tostring(state_names[i]), {state = states[i], pos = {}}, Q[i]))
end

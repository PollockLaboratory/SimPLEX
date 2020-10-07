--[[

   DYNAMIC DISCRETE PARAMETERS.

]]
Model.set_name("JC69 - general")

model_states = Config.get_string_array("MODEL.states")
States.set(model_states)

-- Setup Categories

lower_bound = 0.0
upper_bound = 0.003
step_size = 0.0004

categories = {}

for i = 0, (upper_bound - lower_bound)/step_size do
   l = lower_bound + (i * step_size)
   w = lower_bound + ((i + 1) * step_size)
   init = (l + w)/2.0

   categories[i+1] = Parameter.new("Cat-" .. i, "continuous", {initial_value = init, step_size = 0.01, lower_bound = 0})
end

-- Set lower bounds.
for i = 2, #categories do
   categories[i]:set_lower_bound(categories[i-1])
end

-- Set upper bounds.
for i = 1, #categories-1 do
   categories[i]:set_upper_bound(categories[i+1])
end

rate_categories = Categories.new("RateCategories", categories)

x = Parameter.new("x", "discrete", {categories = rate_categories})

for i=1,#model_states do
	rv_list = {}
	for j=1,#model_states do
	   if i ~= j then
	      rv_list[j] = x
	   else
	      rv_list[j] = Parameter.new("virtual-"..tostring(States[i]), "virtual", {})
	   end
	end
	Model.add_rate_vector(RateVector.new("RV-"..tostring(States[i]), {state = States[i], pos = {}}, rv_list))
end

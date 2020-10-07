--[[

   DISCRETE PARAMETERS.

]]

Model.set_name("JC69 - general")

-- Set the states.
model_states = Config.get_string_array("MODEL.states")
States.set(model_states)

-- Create a table of rate categories given options in the config file.
n = Config.get_float("MODEL.number_of_categories")

lower_limit = Config.get_float("MODEL.lower_limit")
upper_limit = Config.get_float("MODEL.upper_limit")

step_size = (upper_limit - lower_limit)/(n - 1)

categories = {}
for i=0,n-1 do
   table.insert(categories, lower_limit + (i * step_size))
end

-- Create Categories object.
rate_categories = Categories.new("RateCats", categories)

-- Create the discrete parameter.
x = Parameter.new("x", "discrete", {categories = rate_categories})

-- Create all the rate vectors for the substitution model.
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

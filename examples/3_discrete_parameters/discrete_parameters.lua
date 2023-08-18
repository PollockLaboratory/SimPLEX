--[[

   DISCRETE PARAMETERS.
   This model is similar to the structure of the previous model except rather than treating the rate
   parameters as continuous the parameter is sampled from a set of predetermined categories This has
   a similar effect to reducing the degrees of freedom in the substitution model, as the size of the
   state space the Markov chain has to explore is reduced. In SimPLEX parameters that are sampled from
   a finite set of categories are referred to a "discrete".

   In this model the rate categories themselves are not sampled.
]]

Model.set_name("JC69 - general")

-- Set the states.
model_states = Config.get_string_array("MODEL.states")

States.new("nucleotide", model_states,
	   {sequences_output = Config.get_str("MODEL.sequences_out_file"), substitutions_output = Config.get_str("MODEL.substitutions_out_file")})

Data.load_state("nucleotide", Config.get_str("MODEL.sequences_file"))

-- To start we must construct a table that represents the possible values a parameter can be.
-- We start by reading some configuration parameters from the options file, including the number
-- of categories, n and the range of the categories.
n = Config.get_float("MODEL.number_of_categories")

lower_limit = Config.get_float("MODEL.lower_limit")
upper_limit = Config.get_float("MODEL.upper_limit")

step_size = (upper_limit - lower_limit)/(n - 1)

-- We then create the actual lua table and fill it with the possible parameter values.
categories = {}
for i=0,n-1 do
   table.insert(categories, lower_limit + (i * step_size))
end

-- We then load these categories into SimPLEX using the following:
rate_categories = Categories.new("RateCats", categories)

-- We now can create parameters that are sampled from this table we have just create using the
-- Parameter.new function used in the previous example. However, in this case we set the sampling
-- type to "discrete" (the second argument).
x = Parameter.new("x", "discrete", {categories = rate_categories})

-- Finally, we can also construct the rate vectors in the same way as the previous example.
for i=1,#model_states do
	rv_list = {}
	for j=1,#model_states do
	   if i ~= j then
	      rv_list[j] = x
	   else
	      rv_list[j] = Parameter.new("virtual-"..tostring(States.nucleotide[i]), "virtual", {})
	   end
	end
	Model.add_rate_vector(RateVector.new("RV-"..tostring(States.nucleotide[i]), {domain = "nucleotide", state = States.nucleotide[i], pos = {}}, rv_list))
end

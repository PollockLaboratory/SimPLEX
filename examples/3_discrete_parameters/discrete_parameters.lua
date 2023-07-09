--[[

   DISCRETE PARAMETERS.
   This model is similar to the previous model except rather than treating the rate parameters
   as continuous the parameters are treated as discrete. This has the effect of constraining the
   size of the state space the Markov chain has to explore.

   In this model the rate categories themselves are not sampled.
]]

Model.set_name("JC69 - general")

-- Set the states.
model_states = Config.get_string_array("MODEL.states")

States.new("nucleotide", model_states,
	   {sequences_output = Config.get_str("MODEL.sequences_out_file"), substitutions_output = Config.get_str("MODEL.substitutions_out_file")})

Data.load_state("nucleotide", Config.get_str("MODEL.sequences_file"))

-- Create a lua table of rate categories given options in the config file.
-- This creates a list of values from which the rate parameters will be sampled.
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
	      rv_list[j] = Parameter.new("virtual-"..tostring(States.nucleotide[i]), "virtual", {})
	   end
	end
	Model.add_rate_vector(RateVector.new("RV-"..tostring(States.nucleotide[i]), {domain = "nucleotide", state = States.nucleotide[i], pos = {}}, rv_list))
end

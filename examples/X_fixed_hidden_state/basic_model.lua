--[[

   A VERY SIMPLE SUBSTITUTION MODEL - for nucleotide sequences.
 
]]

-- Set the name of the model - this has not effect on the structure of the model.
Model.set_name("JC69 - nucleotide")

-- Set the state of the model.
model_states = {"a", "t", "c", "g"}
States.new("nucleotide", model_states, {})

Data.load_state("nucleotide", Config.get_str("MODEL.sequences_file"))

--States.new("secondary", { "A", "B" }, {})
--Data.load_site_static_state("secondary", "")

x = Parameter.new("x", "continuous", {initial_value = 0.001, step_size = Config.get_float("MODEL.step_size"), lower_bound = 0.0})

-- Note: the step_size parameter is set by looking in the options TOML file, enabling model paramters to be tweaked without changing
-- the lua script.

-- Before creating the first rate vector, RV_a, we must first create the necessary virtual rate parameter. v_a.
-- To create the virtual substitution rate we use the Parameter.new function again, but this time with the type 'virtual'.
v_a = Parameter.new("virtual-a", "virtual", {})

-- Next we create the substitution rate vector itself using the RateVector.new function. This function has the
-- following three arguments:
-- - name - STRING - specifies the name of the rate vector.
-- - application - TABLE - specifies what conditions the substitution rate vector applies to. Most importantly, this specifies
--                         which state the rate vector applies to.
-- - rate parameters - TABLE OF RATE PARAMETERS - specifies the substituion-rate parameters that make up the rate vector.
RV_a = RateVector.new("RV-a", {domain = "nucleotide", state = "a", pos = {}}, {v_a, x, x, x})

-- Using the same process as to create the first substitution-rate vector, RV_a, we can create the final three substitution-rate vectors.
v_t = Parameter.new("virtual-t", "virtual", {})
RV_t = RateVector.new("RV-t", {domain = "nucleotide", state = "t", pos = {}}, {x, v_t, x, x})

v_c = Parameter.new("virtual-c", "virtual", {})
RV_c = RateVector.new("RV-c", {domain = "nucleotide", state = "c", pos = {}}, {x, x, v_c, x})

v_g = Parameter.new("virtual-g", "virtual", {})
RV_g = RateVector.new("RV-g", {domain = "nucleotide", state = "g", pos = {}}, {x, x, x, v_g})

-- Finally, each of the rate vectors must be explicitally added to the complete substitution model using the model.add_rate_vector function.
Model.add_rate_vector(RV_a)
Model.add_rate_vector(RV_t)
Model.add_rate_vector(RV_c)
Model.add_rate_vector(RV_g)

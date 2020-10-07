--[[

   A VERY SIMPLE SUBSTITUTION MODEL - for nucleotide sequences.

   This is a basic example lua script for a simple substitution model for use with simPLEX.
   This demonstrates all the basic features of simPLEX, including:
   - Setting the states of the model.
   - Initiating parameters.
   - Constructing x vectors, specifying the substitution rate from state i to all other states.
   - Passing those x vectors to simPLEX.

   In this simple example will construct a JC69 (Jukes and Cantor 1969) nucleotide model. This model has a single free
   substitution rate parameter, denoted here as x, and 4 possible state: 'a', 't', 'c', 'g'. This model is expressed as
   the following substitution-rate matrix:

   Q = {q_ij} = {{-3x, x, x, x}, 
                 {x, -3x, x, x},
                 {x, x, -3x, x},
                 {x, x, x, -3x}}

   The diagonal elements of the rate matrix, corresponding to the rate of a state substituting to itself, q_ii, must be
   specificed as a unique virtual substitution rate. Thus the above matrix becomes:

   Q = {q_ij} = {{v_a , x, x, x}, 
                 {x, v_t, x, x},
                 {x, x, v_c, x},
                 {x, x, x, v_g}}

   simPLEX specifies models sets of substitution-rate vectors, which in this case corresponds to each row of the above matrix.
   Thus, four substitution-rate vectors are required for each possible state:
   - RV-a = {v_a, x, x, x}
   - RV-t = {x, v_t, x, x}
   - RV-c = {x, x, v_c, x}
   - RV-g = {x, x, x, v_g}

   Let's now implement this model. 
]]

-- Set the name of the model - this has not effect on the structure of the model.
Model.set_name("JC69 - nucleotide")

-- Set the state of the model.
model_states = {"a", "t", "c", "g"}
States.set(model_states)

-- To create the single free parameter, the Parameter.new function is called, this has three arguments:
--   - name - STRING - specifies the name of the parameter, which will be used in the output files.
--   - type - STRING - specifies the type of the parameter, either 'continuous', discrete', 'fixed' or 'virtual'.
--   - options - TABLE - specifies the unique options for each parameter type.
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
RV_a = RateVector.new("RV-a", {state = "a", pos = {}}, {v_a, x, x, x})

-- Using the same process as to create the first substitution-rate vector, RV_a, we can create the final three substitution-rate vectors.
v_t = Parameter.new("virtual-t", "virtual", {})
RV_t = RateVector.new("RV-t", {state = "t", pos = {}}, {x, v_t, x, x})

v_c = Parameter.new("virtual-c", "virtual", {})
RV_c = RateVector.new("RV-c", {state = "c", pos = {}}, {x, x, v_c, x})

v_g = Parameter.new("virtual-g", "virtual", {})
RV_g = RateVector.new("RV-g", {state = "g", pos = {}}, {x, x, x, v_g})

-- Finally, each of the rate vectors must be explicitally added to the complete substitution model using the model.add_rate_vector function.
Model.add_rate_vector(RV_a)
Model.add_rate_vector(RV_t)
Model.add_rate_vector(RV_c)
Model.add_rate_vector(RV_g)

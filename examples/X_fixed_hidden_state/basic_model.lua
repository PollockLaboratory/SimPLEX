--[[

   A VERY SIMPLE SUBSTITUTION MODEL - for nucleotide sequences.
 
]]

-- Set the name of the model - this has not effect on the structure of the model.
Model.set_name("JC69 - nucleotide")

-- Set the state of the model.
model_states = {"a", "t", "c", "g"}
States.new("nucleotide", model_states, {})

Data.load_state("nucleotide", Config.get_str("MODEL.sequences_file"))

States.new("secondary", { "A", "B" }, {})
Data.load_site_static_state("secondary", "datasets/secondary.fasta")

--x = Parameter.new("x", "continuous", {initial_value = 0.01, step_size = Config.get_float("MODEL.step_size"), lower_bound = 0.0})
x = Parameter.new("x", "fixed", { value = 0.01 })

v_a = Parameter.new("virtual-a", "virtual", {})

RV_a = RateVector.new("RV-a", {domain = "nucleotide", state = "a", secondary="*"}, {v_a, x, x, x})
--v_a2 = Parameter.new("virtual-a2", "virtual", {})
--RV_a2 = RateVector.new("RV-a2", {domain = "nucleotide", state = "a", secondary="B"}, {v_a2, x, x, x})

v_t = Parameter.new("virtual-t", "virtual", {})
RV_t = RateVector.new("RV-t", {domain = "nucleotide", state = "t", secondary="*"}, {x, v_t, x, x})

v_c = Parameter.new("virtual-c", "virtual", {})
RV_c = RateVector.new("RV-c", {domain = "nucleotide", state = "c", secondary="*"}, {x, x, v_c, x})

v_g = Parameter.new("virtual-g", "virtual", {})
RV_g = RateVector.new("RV-g", {domain = "nucleotide", state = "g", secondary="*"}, {x, x, x, v_g})

Model.add_rate_vector(RV_a)
--Model.add_rate_vector(RV_a2)
Model.add_rate_vector(RV_t)
Model.add_rate_vector(RV_c)
Model.add_rate_vector(RV_g)

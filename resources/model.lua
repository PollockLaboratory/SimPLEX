model.set_name("Mixture Model")

states.set(config.get_string_array("MODEL.states"))

param_template = {type = "float", initial_value = 0.01}

NtoN = Parameter.new("NtoN", param_template)
NtoF = Parameter.new("NtoF", param_template)
FtoN = Parameter.new("FtoN", param_template)
FtoF = Parameter.new("FtoF", param_template)

model.add_rate_vector(RateVector.new("Non-functional", {state = states[1], pos = {}}, {NtoN, NtoF}))
model.add_rate_vector(RateVector.new("Functional", {state = states[2], pos = {}}, {FtoN, FtoF}))


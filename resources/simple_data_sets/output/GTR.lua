model.set_name("Mixture Model")

states.set(config.get_string_array("MODEL.states"))

N_freq = Parameter.new("N_freq", "continuous", {initial_value = 0.2, step_size = 0.01, lower_bound = 0.0, upper_bound = 1.0 })

F_freq = Parameter.new("F_freq", "continuous", {initial_value = 0.2, step_size = 0.01, lower_bound = 0.0, upper_bound = 1.0 })

--rate = Parameter.new("rate", "continuous", {initial_value = 0.001, step_size = 0.01, lower_bound = 0.0})

one = Parameter.new("fixed_one", "fixed", {value = 0.01})

NtoN = Parameter.new("NN", "virtual", {})
NtoF = Parameter.named_add("NtoF", N_freq, one)
FtoN = Parameter.named_add("FtoN", F_freq, one)
FtoF = Parameter.new("FF", "virtual", {})

print("Name:", NtoN:type())

model.add_rate_vector(RateVector.new("Non-functional", {state = states[1], pos = {}}, {NtoN, NtoF}))
model.add_rate_vector(RateVector.new("Functional", {state = states[2], pos = {}}, {FtoN, FtoF}))


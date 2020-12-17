Model.set_name("Mixture Model")

States.set(Config.get_string_array("MODEL.states"))

N_freq = Parameter.new("N_freq", "continuous", {initial_value = 0.2, step_size = 0.01, lower_bound = 0.0, upper_bound = 1.0 })

F_freq = Parameter.new("F_freq", "continuous", {initial_value = 0.2, step_size = 0.01, lower_bound = 0.0, upper_bound = 1.0 })

--rate = Parameter.new("rate", "continuous", {initial_value = 0.001, step_size = 0.01, lower_bound = 0.0})

one = Parameter.new("fixed_one", "fixed", {value = 0.01})

NtoN = Parameter.new("NN", "virtual", {})
NtoF = Parameter.named_add("NtoF", N_freq, one)
FtoN = Parameter.named_add("FtoN", F_freq, one)
FtoF = Parameter.new("FF", "virtual", {})

print("Name:", NtoN:type())

Model.add_rate_vector(RateVector.new("Non-functional", {state = States[1], pos = {}}, {NtoN, NtoF}))
Model.add_rate_vector(RateVector.new("Functional", {state = States[2], pos = {}}, {FtoN, FtoF}))


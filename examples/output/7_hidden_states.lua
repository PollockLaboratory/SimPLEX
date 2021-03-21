--[[

   Hidden state.

]]--

Model.set_name("Hidden State.")

-- Base
model_states = Config.get_string_array("MODEL.states")
States.new("primary", model_states, {sequences_output = Config.get_str("MODEL.sequences_out_file"), substitutions_output = Config.get_str("MODEL.substitutions_out_file")})
Data.load_state("primary", Config.get_str("MODEL.sequences_file"))

od_states = {"O", "D"}
States.new("orderVdisorder", od_states, {sequences_output = "od_sequences.fasta", substitutions_output = "od_substitutions.out"})

Data.load_state("orderVdisorder", "data_sets/hidden_state.efasta");

-- Create parameters for equilibrium frequencies.

order_rate = Parameter.new("base-order", "continuous", {initial_value = 0.001, step_size = Config.get_float("MODEL.step_size"), lower_bound = 0.0});
disorder_rate = Parameter.new("base-disorder", "continuous", {initial_value = 0.001, step_size = Config.get_float("MODEL.step_size"), lower_bound = 0.0});

-- equil_freq = {}
Q_order = {}
Q_disorder = {}

for i=1,#model_states do
   Q_order[i] = {}
   Q_disorder[i] = {}
end

for i=1,#model_states do
	for j=i,#model_states do
	   if i ~= j then
	      Q_order[i][j] = order_rate
	      Q_order[j][i] = order_rate

	      Q_disorder[i][j] = disorder_rate
	      Q_disorder[j][i] = disorder_rate
	   else
	      Q_order[i][j] = Parameter.new("order-virtual-"..tostring(States.primary[i]), "virtual", {})
	      Q_disorder[i][j] = Parameter.new("disorder-virtual-"..tostring(States.primary[i]), "virtual", {})
	   end
	end
	Model.add_rate_vector(RateVector.new("RV-order-"..tostring(States.primary[i]), {state = States.primary[i], orderVdisorder="O", pos = {}}, Q_order[i]))
	Model.add_rate_vector(RateVector.new("RV-disorder-"..tostring(States.primary[i]), {state = States.primary[i], orderVdisorder="D", pos = {}}, Q_disorder[i]))
end

OtoD = Parameter.new("OtoD", "continuous", {initial_value = 0.001, step_size = Config.get_float("MODEL.step_size"), lower_bound = 0.0 })

DtoO = Parameter.new("DtoO", "continuous", {initial_value = 0.001, step_size = Config.get_float("MODEL.step_size"), lower_bound = 0.0 })

Model.add_rate_vector(RateVector.new("RV-O",
				     {domain = "orderVdisorder", state = "O", primary="*", pos = {}},
				     {Parameter.new("virtual-O", "virtual", {}), OtoD}))

Model.add_rate_vector(RateVector.new("RV-D",
				     {domain = "orderVdisorder", state = "D", primary="*", pos = {}},
				     {DtoO, Parameter.new("virtual-D", "virtual", {})}))

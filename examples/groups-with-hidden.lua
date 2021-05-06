--[[

   (ANOTHER) VERY SIMPLE SUBSTITUTION MODEL - for primary sequences.

   This script builds on the model presented in the previous example, but rather than using a predefined
   list of states, expands the model to any list of n states defined in the configuration file.

   Thus the new general model is:
   q_ij = x, where i != j
   q_ij = -(n-1)x = v_i, where i == j, where v_i is the virtual substitution rate.

]]

-- Group residues.

Model.set_name("Groups with hidden")

-- Base.
model_states = Config.get_string_array("MODEL.states")
States.new("primary", model_states,
	   {sequences_output = Config.get_str("MODEL.sequences_out_file"), substitutions_output = Config.get_str("MODEL.substitutions_out_file")})

Data.load_state("primary", Config.get_str("MODEL.sequences_file"))

od_states = Config.get_string_array("MODEL.od_states")
States.new("orderVdisorder", od_states, {sequences_output = Config.get_str("MODEL.od_sequences_out_file"), substitutions_output = Config.get_str("MODEL.od_substitutions_out_file")})

Data.load_state("orderVdisorder", Config.get_str("MODEL.od_sequences_file"));

-- Read file to group amino acid residues into groups.
function lines_from(file)
   lines = {}
   for line in io.lines(file) do
      lines[#lines + 1] = line
      print(line)
   end
   return lines
end

groups = lines_from(Config.get_root_directory()..Config.get_str("MODEL.groups_file"))

-- Create rate parameters for each group.
subs_to_rate_order = {}
subs_to_rate_disorder = {}
for i, str in pairs(groups) do
   rate_order = Parameter.new("Group"..tostring(i).."-order", "continuous",
			      {initial_value = 0.001, step_size = Config.get_float("MODEL.step_size"), lower_bound = 0.0})
   rate_disorder = Parameter.new("Group"..tostring(i).."-disorder", "continuous",
				 {initial_value = 0.001, step_size = Config.get_float("MODEL.step_size"), lower_bound = 0.0})
   
   for sub_pair in str:gmatch("%w+") do
      print(i, sub_pair)
      subs_to_rate_order[sub_pair] = rate_order
      subs_to_rate_disorder[sub_pair] = rate_disorder
   end
end

-- Build rate vectors.
for i=1,#model_states do
   rv_list_order = {}
   rv_list_disorder = {}
   for j=1,#model_states do
      if i ~= j then
	 sub_pair = States.primary[i]..States.primary[j]
	 print(sub_pair)
	 rv_list_order[j] = subs_to_rate_order[sub_pair]
	 rv_list_disorder[j] = subs_to_rate_disorder[sub_pair]
      else
	 rv_list_order[j] = Parameter.new("virtual-"..tostring(States.primary[i]).."-order", "virtual", {})
	 rv_list_disorder[j] = Parameter.new("virtual-"..tostring(States.primary[i]).."-disorder", "virtual", {})
      end
   end
   Model.add_rate_vector(RateVector.new("RV-"..tostring(States.primary[i]), {domain = "primary", state = States.primary[i], orderVdisorder = "O", pos = {}}, rv_list_order))
   Model.add_rate_vector(RateVector.new("RV-"..tostring(States.primary[i]), {domain = "primary", state = States.primary[i], orderVdisorder = "D", pos = {}}, rv_list_disorder))
end

print("OD:", Config.get_bool("MODEL.equal_od_rates"));
if Config.get_bool("MODEL.equal_od_rates") == 1 then
   print("Order -> disorder and disorder -> order rates constrained to be equal.")
   OtoD = Parameter.new("OD_transition_rate", "continuous", {initial_value = 0.001, step_size = Config.get_float("MODEL.step_size"), lower_bound = 0.0 })
   DtoO = OtoD
else
   print("Order -> disorder and disorder -> order rates are separate rate parameters.")
   OtoD = Parameter.new("OtoD", "continuous", {initial_value = 0.001, step_size = Config.get_float("MODEL.step_size"), lower_bound = 0.0 })

   DtoO = Parameter.new("DtoO", "continuous", {initial_value = 0.001, step_size = Config.get_float("MODEL.step_size"), lower_bound = 0.0 })
end

Model.add_rate_vector(RateVector.new("RV-O",
				     {domain = "orderVdisorder", state = "O", primary="*", pos = {}},
				     {Parameter.new("virtual-O", "virtual", {}), OtoD}))

Model.add_rate_vector(RateVector.new("RV-D",
				     {domain = "orderVdisorder", state = "D", primary="*", pos = {}},
				     {DtoO, Parameter.new("virtual-D", "virtual", {})}))



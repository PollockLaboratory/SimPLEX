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
model_states = { 'a', 't', 'c', 'g' }
States.new("primary", model_states,
           {sequences_output = Config.get_str("MODEL.sequences_out_file"),
            substitutions_output = Config.get_str("MODEL.substitutions_out_file")})

Data.load_state("primary", Config.get_str("MODEL.sequences_file"))

hidden_states = { 'A', 'B' }
States.new("hidden", hidden_states, {sequences_output = Config.get_str("MODEL.hidden_sequences_out_file"), substitutions_output = Config.get_str("MODEL.hidden_substitutions_out_file")})

Data.load_state("hidden", Config.get_str("MODEL.hidden_sequences_file"));

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
subs_to_rate_A = {}
subs_to_rate_B = {}
for i, str in pairs(groups) do
   rate_A = Parameter.new("Group"..tostring(i).."-A", "continuous",
			      {initial_value = 0.001, step_size = Config.get_float("MODEL.step_size"), lower_bound = 0.0})
   rate_B = Parameter.new("Group"..tostring(i).."-B", "continuous",
				 {initial_value = 0.001, step_size = Config.get_float("MODEL.step_size"), lower_bound = 0.0})
   
   for sub_pair in str:gmatch("%w+") do
      print(i, sub_pair)
      subs_to_rate_A[sub_pair] = rate_A
      subs_to_rate_B[sub_pair] = rate_B
   end
end

-- Build rate vectors.
for i=1,#model_states do
   rv_list_A = {}
   rv_list_B = {}
   for j=1,#model_states do
      if i ~= j then
	 sub_pair = States.primary[i]..States.primary[j]
	 print(sub_pair)
	 rv_list_A[j] = subs_to_rate_A[sub_pair]
	 rv_list_B[j] = subs_to_rate_B[sub_pair]
      else
	 rv_list_A[j] = Parameter.new("virtual-"..tostring(States.primary[i]).."-hiddenA", "virtual", {})
	 rv_list_B[j] = Parameter.new("virtual-"..tostring(States.primary[i]).."-hiddenB", "virtual", {})
      end
   end
   Model.add_rate_vector(RateVector.new("RV-"..tostring(States.primary[i]), {domain = "primary", state = States.primary[i], hidden = "A", pos = {}}, rv_list_A))
   Model.add_rate_vector(RateVector.new("RV-"..tostring(States.primary[i]), {domain = "primary", state = States.primary[i], hidden = "B", pos = {}}, rv_list_B))
end

AtoB = Parameter.new("AtoB", "continuous", {initial_value = 0.001, step_size = Config.get_float("MODEL.step_size"), lower_bound = 0.0 })

BtoA = Parameter.new("BtoA", "continuous", {initial_value = 0.001, step_size = Config.get_float("MODEL.step_size"), lower_bound = 0.0 })

--
Model.add_rate_vector(RateVector.new("RV-A",
				     {domain = "hidden", state = "A", primary="*", pos = {}},
				     {Parameter.new("virtual-A", "virtual", {}), AtoB}))

Model.add_rate_vector(RateVector.new("RV-B",
				     {domain = "hidden", state = "B", primary="*", pos = {}},
				     {BtoA, Parameter.new("virtual-A", "virtual", {})}))



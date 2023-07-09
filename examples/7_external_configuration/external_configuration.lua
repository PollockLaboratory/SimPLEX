--[[

   (ANOTHER) VERY SIMPLE SUBSTITUTION MODEL - for nucleotide sequences.

   This script builds on the model presented in the previous example, but rather than using a predefined
   list of states, expands the model to any list of n states defined in the configuration file.

   Thus the new general model is:
   q_ij = x, where i != j
   q_ij = -(n-1)x = v_i, where i == j, where v_i is the virtual substitution rate.

]]

-- Group residues.

Model.set_name("external_configuration")

-- Set the states.
model_states = Config.get_string_array("MODEL.states")
States.new("nucleotide", model_states,
	   {sequences_output = Config.get_str("MODEL.sequences_out_file"), substitutions_output = Config.get_str("MODEL.substitutions_out_file")})

Data.load_state("nucleotide", Config.get_str("MODEL.sequences_file"))

-- Read file to group amino acid residues into groups.
function lines_from(file)
   lines = {}
   for line in io.lines(file) do
      lines[#lines + 1] = line
      print(line)
   end
   return lines
end

groups = lines_from(Config.get_root_directory().."/datasets/nucleotide_groups.txt")

-- Create rate parameters for each group.
subs_to_rate = {}
for i, str in pairs(groups) do
   rate = Parameter.new("Group"..tostring(i), "continuous",
			{initial_value = 0.001, step_size = Config.get_float("MODEL.step_size"), lower_bound = 0.0})
   for sub_pair in str:gmatch("%w+") do
      print(i, sub_pair)
      subs_to_rate[sub_pair] = rate
   end
end

-- Build rate vectors.
for i=1,#model_states do
	rv_list = {}
	for j=1,#model_states do
	   if i ~= j then
	      sub_pair = States.nucleotide[i]..States.nucleotide[j]
	      print(sub_pair)
	      rv_list[j] = subs_to_rate[sub_pair]
	   else
	      rv_list[j] = Parameter.new("virtual-"..tostring(States.nucleotide[i]), "virtual", {})
	   end
	end
	Model.add_rate_vector(RateVector.new("RV-"..tostring(States.nucleotide[i]), {domain = "nucleotide", state = States.nucleotide[i], pos = {}}, rv_list))
end


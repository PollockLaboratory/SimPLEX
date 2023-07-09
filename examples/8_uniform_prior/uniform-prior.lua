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

n = Config.get_float("MODEL.n_hidden_categories")
hidden_states = {}
for i=1,n do
   table.insert(hidden_states, string.char(64+i))
end

States.new("hidden_category", hidden_states, {sequences_output = Config.get_str("MODEL.od_sequences_out_file"), substitutions_output = Config.get_str("MODEL.od_substitutions_out_file")})

Data.create_uniform_prior("hidden_category");

-- Read file to group amino acid residues into groups.
function lines_from(file)
   lines = {}
   for line in io.lines(file) do
      lines[#lines + 1] = line
      print(line)
   end
   return lines
end

param_template = {initial_value = 0.001, step_size = Config.get_float("MODEL.step_size"), lower_bound = 0.0}
-- Build rate vectors.
for k, hidden_state in pairs(hidden_states) do
   for i=1,#model_states do
      rv_list = {}
      for j=1,#model_states do
	 if i ~= j then
	    sub_pair = States.primary[i]..States.primary[j]
	    print(sub_pair)
	    rv_list[j] = Parameter.new(sub_pair.."-"..hidden_state, "continuous", param_template)
	 else
	    rv_list[j] = Parameter.new("virtual-"..tostring(States.primary[i]).."-"..hidden_state, "virtual", {})
	 end
      end
      Model.add_rate_vector(RateVector.new("RV-"..tostring(States.primary[i]), {domain = "primary", state = States.primary[i], hidden_category = hidden_state, pos = {}}, rv_list))
   end
end

-- HIDDEN STATE RATE VECTORS.
function tail(i, list)
   return { select(i+1, table.unpack(list)) }
end

hidden_transition_pairs = {}
for i, s1 in pairs(hidden_states) do
   for j, s2 in pairs(tail(i, hidden_states)) do
      hidden_transition_rate = Parameter.new(s1.."to"..s2, "continuous",  {initial_value = 0.001, step_size = Config.get_float("MODEL.step_size"), lower_bound = 0.0 })
      hidden_transition_pairs[s1..s2] = hidden_transition_rate
      hidden_transition_pairs[s2..s1] = hidden_transition_rate
      print(s1..s2, s2..s1)
   end
end

for i=1,#hidden_states do
   rv_list = {}
   print(i, j)
   for j=1,#hidden_states do
      if i ~= j then
	 sub_pair = States.hidden_category[i]..States.hidden_category[j]
	 rv_list[j] = hidden_transition_pairs[sub_pair]
      else
	 rv_list[j] = Parameter.new("virtual-"..States.hidden_category[i], "virtual", {})
      end 
   end
   Model.add_rate_vector(RateVector.new("RV-"..States.hidden_category[i],
					{domain = "hidden_category", state = States.hidden_category[i], primary="*", pos = {}},
					rv_list))
end

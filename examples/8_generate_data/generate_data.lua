--[[

   This example shows how a uniformly distributed data can be generate for the hidden states at the tree tips, rather
   than providing input data explicitly. We demonstrate this here with substitution model that has N hidden states,
   however, we have no data about which of these hidden states the residues in the tip sequences might be. Thus, we set
   the probability of each site being in each of the hidden states to 1/N, using the Data.generate_uniform_data function.
]]

-- Nucleotide state.
model_states = {'a', 't', 'c', 'g'}
States.new("nucleotide", model_states,
	   {sequences_output = Config.get_str("MODEL.sequences_out_file"), substitutions_output = Config.get_str("MODEL.substitutions_out_file")})

Data.load_state("nucleotide", Config.get_str("MODEL.sequences_file"))

n = Config.get_float("MODEL.n_hidden_categories")
hidden_states = {}
for i=1,n do
   table.insert(hidden_states, string.char(64+i))
end

States.new("hidden", hidden_states, {sequences_output = Config.get_str("MODEL.od_sequences_out_file"), substitutions_output = Config.get_str("MODEL.od_substitutions_out_file")})

-- Here, rather than loading or reading from a source file, the data is generated to match the pattern (in terms of gaps) of the previously
-- defined state - in this case the nucleotide state. To reiterate each residue in a sequence is given an equal probability of being in each hidden state.
Data.generate_uniform_data("hidden");

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
	    sub_pair = States.nucleotide[i]..States.nucleotide[j]
	    print(sub_pair)
	    rv_list[j] = Parameter.new(sub_pair.."-"..hidden_state, "continuous", param_template)
	 else
	    rv_list[j] = Parameter.new("virtual-"..tostring(States.nucleotide[i]).."-"..hidden_state, "virtual", {})
	 end
      end
      Model.add_rate_vector(RateVector.new("RV-"..tostring(States.nucleotide[i]), {domain = "nucleotide", state = States.nucleotide[i], hidden = hidden_state, pos = {}}, rv_list))
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
         sub_pair = States.hidden[i]..States.hidden[j]
         rv_list[j] = hidden_transition_pairs[sub_pair]
      else
         rv_list[j] = Parameter.new("virtual-"..States.hidden[i], "virtual", {})
      end 
   end
   Model.add_rate_vector(RateVector.new("RV-"..States.hidden[i],
					{domain = "hidden", state = States.hidden[i], nucleotide="*", pos = {}},
					rv_list))
end

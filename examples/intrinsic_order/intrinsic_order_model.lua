--[[
   
   This is the final model used in the paper:
   "Disentangling multi-dimensional context-dependent susbtitution patterns that vary through time." by Pike et al

   This substitution model aims to distinguish the difference in amino acid substitutions patterns between
   intrinsic order and intrinsic disorder. For each intrinsic order context, there are separate non-reversible amino acid substitution
   rate matrices. Naively, this would result in a model containing 760 free parameters - an excessive amount given many of the amino acid
   pairs do not provide any useful information to descriminate intrinsic order contexts. From previous analyses we have
   determined amino acid pair that either do not occur in our dataset, and are lablled as rare, or provide little descriminatory
   information, labelled as uninformative. We use this to simplify the model and merge substitution rate parameters.
   
   See paper for more details.

]]

-- Utility functions to assist reading in text files containing lists of amino acid pairs.
function file_exists(file)
  local f = io.open(file, "rb")
  if f then f:close() end
  return f ~= nil
end

-- get all lines from a file, returns an empty 
-- list/table if the file does not exist
function lines_from(file)
   if not file_exists(file) then
      print("NOT EXIST")
      return {}
   end
   local lines = {}
   for line in io.lines(file) do 
      lines[#lines + 1] = line
   end
   return lines
end

-- CONFIGURE STATE DIMENSIONS
-- There are two state dimensions in this model: "amino" representing amino acid substitutions and "hidden_intrinsic_order" representing
-- substitutions between intrinsic order context. Dimensions are referred to a domains in the SimPLEX API.

-- AMINO ACID DIMENSION
model_states = Config.get_string_array("MODEL.states")
States.new("amino", model_states,
           {sequences_output = Config.get_str("MODEL.sequences_out_file"), substitutions_output = Config.get_str("MODEL.substitutions_out_file")})

Data.load_state("amino", Config.get_str("MODEL.sequences_file"))

-- INTRINSIC ORDER DIMENSION
intrinsic_order_states = Config.get_string_array("MODEL.intrinsic_order_states")
States.new("hidden_intrinsic_order", intrinsic_order_states, {sequences_output = Config.get_str("MODEL.intrinsic_order_sequences_out_file"), substitutions_output = Config.get_str("MODEL.intrinsic_order_substitutions_out_file")})

Data.load_state("hidden_intrinsic_order", Config.get_str("MODEL.intrinsic_order_sequences_file"));

-- Load and parse the data files that contain the list amino acids that have been categorized into the uninformative and rare groups.
-- Any amino acids not listed in these two files are inferred to be in the informative category.

-- RARE SUBSTITUTIONS
print("Rare substitutions:")
rare_subs = {}
c = 0
for k, v in pairs(lines_from(Config.get_str("MODEL.rare_substitutions_file"))) do
  io.write(v, " ")
  if k % 15 == 0 then
	  io.write("\n")
	  io.flush()
  end
  c = c + 1
  rare_subs[v] = true
end
io.write("\nTotal rare amino acid substitutions:", c)
io.write("\n\n")
io.flush()

-- UNINFORMATIVE SUBSTITUTIONS
print("Uninformative substitutions:")
uninformative_subs = {}
c = 0
for k, v in pairs(lines_from(Config.get_str("MODEL.uninformative_substitutions_file"))) do
  io.write(v, " ")
  if k % 15 == 0 then
	  io.write("\n")
	  io.flush()
  end
  c = c + 1
  uninformative_subs[v] = true
end
io.write("\nTotal uninformative amino acid substitutions:", c)
io.write("\n\n")
io.flush()


function amino_acid_pair_from_index(i, j)
   return(tostring(States.amino[i])..tostring(States.amino[j]))
end

-- Amino acid pairs that are categorized as rare, have a single fixed slow rate.
-- The non-sampled pseudo-parameter for this fixed rate is created here.
rare_rate = Parameter.new("rare", "fixed", {value = 0.0001})

init_val_order = 0.001
init_val_disorder = 0.001

-- Construct the array of arrays that will become the amino acid rate matrices in either the intrinsically
-- ordered context or the intrinsically disordered context.
Q_order = {}
Q_disorder = {}

for i=1,#model_states do
   Q_order[i] = {}
   Q_disorder[i] = {}
end

-- Counts of each amino acid category.
n_rare = 0
n_uninformative = 0
n_informative = 0

function create_parameters(i, j)
   -- This function takes an amino acid pair specified by their index and creates and returns the necessary rate parameters.
   -- Each amino acid pair is checked to which category it is in, either: Informative, Uninformative, Rare.
   -- This returns a table with two keys: order, disorder where the corresponding vales contain the rate parameters for the specified amino acid.
   local aa_name = amino_acid_pair_from_index(i, j)
   if (rare_subs[aa_name] == true) then
      -- RARE
      -- the previously defined rare rate parameter is returned for both contexts.
      n_rare = n_rare + 1
      return({order = rare_rate, disorder = rare_rate})
   elseif(uninformative_subs[aa_name] == true) then
      -- UNINFORMATIVE
      -- a single new rate parameter is constructed that is shared across both contexts..
      n_uninformative = n_uninformative + 1
      local s = Parameter.new(amino_acid_pair_from_index(i, j).."-shared", "continuous",
                              {initial_value = 0.0001, step_size = Config.get_float("MODEL.step_size"), lower_bound = 0.0})
      return({order = s, disorder = s})
   else
      -- INFORMATIVE
      -- a two new rate parameters are constructed for the intrinsically ordered and intrinsically disordered contexts.
      n_informative = n_informative + 1
      local o = Parameter.new(amino_acid_pair_from_index(i, j).."-order", "continuous",
                              {initial_value = 0.0001, step_size = Config.get_float("MODEL.step_size"), lower_bound = 0.0})
      local d = Parameter.new(amino_acid_pair_from_index(i, j).."-disorder", "continuous",
                              {initial_value = 0.0001, step_size = Config.get_float("MODEL.step_size"), lower_bound = 0.0})
      return({order = o, disorder = d})
   end
end

-- Construct the rate vectors and rate matrices for the amino acid substitutions.
for i=1,#model_states do
   for j=i,#model_states do
      if i ~= j then
         p_ij = create_parameters(i, j)
         Q_order[i][j] = p_ij.order
         Q_disorder[i][j] = p_ij.disorder
	 
         p_ji = create_parameters(j, i)
         Q_order[j][i] = p_ji.order
         Q_disorder[j][i] = p_ji.disorder
      else
         -- Virtual Substitution rate.
         Q_order[i][j] = Parameter.new("order-virtual-"..tostring(States.amino[i]), "virtual", {})
         Q_disorder[i][j] = Parameter.new("disorder-virtual-"..tostring(States.amino[i]), "virtual", {})
      end
   end

   -- Add rate vectors to model.
   Model.add_rate_vector(RateVector.new("RV-order-"..tostring(States.amino[i]),
                                        {domain = "amino", state = States.amino[i], hidden_intrinsic_order="O", pos = {}}, Q_order[i]))
   Model.add_rate_vector(RateVector.new("RV-disorder-"..tostring(States.amino[i]),
                                        {domain = "amino", state = States.amino[i], hidden_intrinsic_order="D", pos = {}}, Q_disorder[i]))
end

io.write("Total by group: Informative: ", n_informative, " Uninformative: ", n_uninformative, " Rare: ", n_rare, "\n")
io.flush()

-- Construct the rate vectors and rate matrices for the intrinsic disorder context.
-- A fixed transition rate is set as the model becomes very unstable otherwise.
intrinsic_order_transition = Parameter.new("intrinsic-order-transition", "fixed", { value = Config.get_float("MODEL.intrinsic_order_transition") })

Model.add_rate_vector(RateVector.new("RV-O",
				     {domain = "hidden_intrinsic_order", state = "O", amino="*", pos = {}},
				     {Parameter.new("virtual-O", "virtual", {}), intrinsic_order_transition}))

Model.add_rate_vector(RateVector.new("RV-D",
				     {domain = "hidden_intrinsic_order", state = "D", amino="*", pos = {}},
				     {intrinsic_order_transition, Parameter.new("virtual-D", "virtual", {})}))

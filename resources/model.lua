model.set_name("Mixture Model")

states.set(config.get_string_array("MODEL.states"))
parameter_type = config.get_str("MODEL.parameter_type")

if parameter_type == "continuous" then
   param_template = {initial_value = 0.05, step_size = 0.01}
elseif parameter_type == "discrete" then
   cats = {}
   val = 0.00
   for i=1,config.get_int("MODEL.num_categories") do
      val = val + 0.01
      cats[i] = val
      print(val)
   end
   print(#cats)
   for k,v in pairs(cats) do
      print(v)
   end 
   param_template = {categories = cats}
end

NtoN = Parameter.new("NtoN", parameter_type, param_template)
NtoF = Parameter.new("NtoF", parameter_type, param_template)
FtoN = Parameter.new("FtoN", parameter_type, param_template)
FtoF = Parameter.new("FtoF", parameter_type, param_template)

model.add_rate_vector(RateVector.new("Non-functional", {state = states[1], pos = {}}, {NtoN, NtoF}))
model.add_rate_vector(RateVector.new("Functional", {state = states[2], pos = {}}, {FtoN, FtoF}))


### Summary
## Biofuel supply chain
# 2 step process
# extraction: seeds -> vegtable oil (calculated with crop_yields AND crop_oil_contents)
# transesterification: vegrable oil + methanol -> biodiesel
#   0.9l biodiesel <- 1l veetable oil, 0.2l methanol
# later: refining: biodiesel + petrol diesel -> refined biodiesel

#
## Crops : soya, sunflower, cutton
# Limited area (1600 ha)
# Different yields(tonnes per ha)
# Different water demands (megalitres per ha (Ml/ha))
# Water limited to 5000 Ml
#
## Final products: B5, B30, B100
# Final: refined biodiesel + petrol diesel -> final fuel
#   Diffirent concentrations for different customers
# Limited amount of petrol diesel: 150_000l
# biodiesel + (maybe petrol diesel) -> standard diesel engines
#   Blends are most commonlyly used for retail fuel marketplace
#   Blends are denoted its B factor which states athe amount of biodiesel
#   Blends: B5, B30, B100
#     B5 "usual diesel"(in eu) is the most common
# Demand: 280_000l
#
### Model
## Optimize: profit from biodiesel production
# profit = sum of each each product: amount sold * (product_price - product_tax) - production cost
## Variables:
# x -- crop area (ha)
# y -- product liters (l)

using JuMP, Clp

arg = length(ARGS) > 0 ? ARGS[1] : ""


### Data:
petrol_diesel_limit = 150_000 # Liters of available petrol disel
## price (Euro per liter)
price_methanol = 1.5 # (Euro/l)
price_petrol_diesel = 1 # (Euro/l)
if arg == "d"
    price_petrol_diesel = 1.2
end

## Crops: soybeans, sunflower seeds, cotton seeds
crop_total_area = 1_600 # ha
crop_total_water_limit = 5_000 # Ml
crop_count = 3
crop_yields = 1000.0 .* [2.6, 1.4, 0.9] # 1000 * t/ha = kg/ha
if arg == "c" # Because the big differentiating factor that makes sunflower seeds less profitable than soybeans is the low yield, other than that the numbers are either similar or better for sunflower seeds than soybeans.
    crop_yields = 1000.0 .* [2.6, 2.1, 0.9]
end
crop_water_demands = [5.0, 4.2, 1.0] # Ml/ha
crop_oil_contents = [0.178, 0.216, 0.433] # l/kg

## Final products: B5, B30, B100
prod_total_demand = 280_000 # l
prod_count = 3
prod_biodiesel_ratios = [0.05, 0.3, 1.0]
prod_petrol_diesel_ratios = 1.0 .- prod_biodiesel_ratios
prod_prices = [1.43, 1.29, 1.16] # Euro/l
prod_taxes = [0.20, 0.05, 0.0]
if arg == "e"
    prod_taxes = [0.12, 0.05, 0.0]
end


### Model:
model = Model(Clp.Optimizer)
@variable(model, x[1:crop_count], lower_bound=0) # ha
@variable(model, y[1:prod_count], lower_bound=0) # l
if arg == "a1"
    petrol_diesel_limit = @variable(model, a_var)
elseif arg == "a2"
    crop_total_water_limit = @variable(model, a_var)
elseif arg == "a3"
    crop_total_area = @variable(model, a_var)
end

water_usage         = sum(x .* crop_water_demands) # ha*(Ml/ha) = Ml 
# Extraction of vegtable oils from seeds
# units: ha * (kg/ha) * (l/kg) = l 
vegtable_oil       = sum(x .* crop_yields .* crop_oil_contents) # l 
# 0.9l biodiesel <- 1l veetable oil, 0.2l methanol
biodiesel_produced  = 0.9 * vegtable_oil  # l
methanol_usage      = 0.2 * vegtable_oil # l
# Each product is a mix of biodiesel and petrol diesel
biodiesel_usage     = sum(y .* prod_biodiesel_ratios) # l
petrol_diesel_usage = sum(y .* prod_petrol_diesel_ratios) # l
@assert all(isapprox.(prod_biodiesel_ratios .+ prod_petrol_diesel_ratios, 1.0)) # They must add up to 1

# Crop constraints
constraint_area = @constraint(model, sum(x) <= crop_total_area) # ha <= ha
constraint_water = @constraint(model, water_usage <= crop_total_water_limit) # Ml <= Ml
# Product constraints
# @constraint(model, sum(y) == prod_total_demand) # l == l
constraint_demand = @constraint(model, sum(y) >= prod_total_demand) # l >= l
constraint_petrol_diesel = @constraint(model, petrol_diesel_usage <= petrol_diesel_limit) # l <= l
# Constraint between x and y: We can only make so much of each product from the biodiesel
constraint_biodiesel = @constraint(model, biodiesel_usage <= biodiesel_produced) # l <= l

cost_methanol = price_methanol * methanol_usage # (Euro/l) * l = Euro
cost_petrol_diesel = price_petrol_diesel * petrol_diesel_usage # (Euro/l) * l = Euro
costs = cost_methanol .+ cost_petrol_diesel # Euro + Euro = Euro

revenues     = y .* prod_prices       # l * (Euro/l) = Euro
taxes        = revenues .* prod_taxes # Euro
profit = sum(revenues) - sum(taxes) - costs # Euro - Euro - Euro = Euro

if arg == "a1" || arg == "a2" || arg == "a3"
    @objective(model, Min, a_var)
else
    @objective(model, Max, profit)
end
println("Model: ", model)
optimize!(model)

println("Termination status: ", termination_status(model))
println("primal_status: ", primal_status(model))
println("objective_value: ", objective_value(model))
println("x: ", value.(x))
println("y: ", value.(y))
println("sum(x): ", sum(value.(x)))
println("sum(y): ", sum(value.(y)))
println("water_usage: ", value(water_usage))
println("vegtable_oil: ", value(vegtable_oil))
println("methanol_usage: ", value(methanol_usage))
println("biodiesel_produced: ", value(biodiesel_produced))
println("biodiesel_usage: ", value(biodiesel_usage))
println("petrol_diesel_usage: ", value(petrol_diesel_usage))
println("cost_methanol: ", value(cost_methanol))
println("cost_petrol_diesel: ", value(cost_petrol_diesel))
println("costs: ", value(costs))
println("revenue: ", value(sum(revenues)))
println("tax: ", value(sum(taxes)))
println("profit: ", value(profit))

if arg == "a1" || arg == "a2" || arg == "a3"
    println("# 3. (a)")
    println("task a variable: ", value(a_var))
end

if arg == "b"
    println("# 3. (b)")
    println("constraint_petrol_diesel: ", constraint_petrol_diesel)
    println("constraint_water: ", constraint_water)
    println("constraint_area: ", constraint_area)
end

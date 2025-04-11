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


### Data:
petrol_diesel_limit = 150_000 # Liters of available petrol disel
## price (Euro per liter)
price_methanol = 1.5 # (Euro/l)
price_petrol_diesel = 1 # (Euro/l)

## Crops: soybeans, sunflower seeds, cotton seeds
crop_total_area = 1_600 # ha
crop_total_water_limit = 5_000 # Ml
crop_count = 3
crop_yields = 1000.0 .* [2.6, 1.4, 0.9] # 1000 * t/ha = kg/ha
crop_water_demands = [5.0, 4.2, 1.0] # Ml/ha
crop_oil_contents = [0.178, 0.216, 0.433] # l/kg

## Final products: B5, B30, B100
prod_total_demand = 280_000 # l
prod_count = 3
prod_biodiesel_ratios = [0.05, 0.3, 1.0]
prod_petrol_diesel_ratios = 1.0 .- prod_biodiesel_ratios
prod_prices = [1.43, 1.29, 1.16] # Euro/l
prod_taxes = [0.20, 0.05, 0]


### Model:
model = Model(Clp.Optimizer)
@variable(model, x[1:crop_count], lower_bound=0) # ha
@variable(model, y[1:prod_count], lower_bound=0) # l

water_usage         = x .* crop_water_demands # ha*(Ml/ha) = Ml 
# Extraction of vegtable oils from seeds
# units: ha * (kg/ha) * (l/kg) = l 
vegtable_oils       = x .* crop_yields .* crop_oil_contents # l 
# 0.9l biodiesel <- 1l veetable oil, 0.2l methanol
biodiesel_produced  = 0.9 .* vegtable_oils  # l
methanol_usage      = 0.2 .* vegtable_oils # l
# Each product is a mix of biodiesel and petrol diesel
biodiesel_usage     = y .* prod_biodiesel_ratios # l
petrol_diesel_usage = y .* prod_petrol_diesel_ratios # l
@assert all(isapprox.(prod_biodiesel_ratios .+ prod_petrol_diesel_ratios, 1.0)) # They must add up to 1

# Crop constraints
@constraint(model, sum(x) <= crop_total_area) # ha <= ha
@constraint(model, sum(water_usage) <= crop_total_water_limit) # Ml <= Ml
# Product constraints
@constraint(model, sum(y) >= prod_total_demand) # l >= l
@constraint(model, sum(petrol_diesel_usage) <= petrol_diesel_limit) # l <= l
# Constraint between x and y: We can only make so much of each product from the biodiesel
@constraint(model, sum(biodiesel_usage) <= sum(biodiesel_produced)) # l <= l

total_cost_methanol = sum(price_methanol .* methanol_usage) # (Euro/l) * l = Euro
total_cost_petrol_diesel = sum(price_petrol_diesel .* petrol_diesel_usage) # (Euro/l) * l = Euro
total_cost = total_cost_methanol + total_cost_petrol_diesel # Euro + Euro = Euro

revenues = y .* prod_prices # l * (Euro/l) = Euro
total_revenue = sum(revenues) # Euro
taxes = revenues .* prod_taxes # Euro
total_tax = sum(taxes) # Euro

@objective(model, Max, total_revenue - total_tax - total_cost) # Euro - Euro - Euro = Euro
optimize!(model)

println("Termination status: ", termination_status(model))
println("primal_status: ", primal_status(model))
println("objective_value: ", objective_value(model))
println("x: ", value.(x))
println("y: ", value.(y))
println("sum(x): ", sum(value.(x)))
println("sum(y): ", sum(value.(y)))

### Summary
## Biofuel supply chain
# 2 step process
# extraction: seeds + ? -> vegtable oil 
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
# x -- crop area???
# y -- product liters???


### Data:

## Biofuels supply chain
# 0.9l unrefined biodisel <- 1l vegtable oil + 0.2l methanol
# 

petrol_diesel_limit = 150_000 # Liters of available petrol disel
## price (Euro per liter)
# methanol, petrol diesel
# fuel_prices = [1.5, 1.0]
price_methanol = 1.5 # (Euro per liter)
price_petrol_diesel = 1 # (Euro per liter)


## Crops
# soybeans, sunflower seeds, cotton seeds
crop_area = 1_600 # ha
crop_water_limit = 5_000 # Mli
crop_count = 3
crop_yields = 1000.0 .* [2.6, 1.4, 0.9] # l/ha
crop_water_demands = [5.0, 4.2, 1.0] # Ml/ha
crop_oil_contents = [0.178, 0.216, 0.433] # l/kg

## Final products: B5, B30, B100
product_demand = 280_000
product_count = 3
product_biodiesel_ratios = [0.05, 0.3, 1.0]
product_petrol_diesel_ratios = 1.0 .- product_biodiesel_ratios
product_prices = [1.43, 1.29, 1.16]
product_taxes = [0.20, 0.05, 0]


using JuMP, Clp

model = Model(Clp.Optimizer)
@variable(model, x[1:crop_count], lower_bound=0)
@variable(model, y[1:product_count], lower_bound=0)

# Crop constraints
@constraint(model, sum(x) <= crop_area)
@constraint(model, sum(x .* crop_water_demands) <= crop_water_limit)

# Product constraints
@constraint(model, sum(y .* product_petrol_diesel_ratios) <= petrol_diesel_limit)
@constraint(model, sum(y) >= product_demand)

# Constraint between x and y: We can only make so much of each product from the biodiesel
@constraint(model, sum(y .* product_biodiesel_ratios) <= sum(x .* crop_yields .* crop_oil_contents .* 0.9))

sold = sum(y .* (product_prices .* (1 .- product_taxes)))
cost_methanol = price_methanol * 0.2 * sum(x .* crop_yields .* crop_oil_contents)
cost_petrol_diesel = price_petrol_diesel .* sum(y .* product_petrol_diesel_ratios) 
cost = cost_methanol + cost_petrol_diesel

@objective(model, Max, sold - cost)
optimize!(model)

println("Termination status: ", termination_status(model))
println("primal_status: ", primal_status(model))
println("objective_value: ", objective_value(model))
println("x: ", value.(x))
println("y: ", value.(y))

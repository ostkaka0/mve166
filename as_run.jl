#!/usr/bin/julia
using JuMP
#using Cbc
using Gurobi
using SparseArrays

arg = length(ARGS) > 0 ? ARGS[1] : ""

if arg == "1a"
    include("as_dat_large.jl")
else
    include("as_dat_small.jl")
end

if arg == "1a"
    T = 125
end

println("# Data:")
println("Components:", Components)
println("T:", T)
println("d:", d)
println("c:", summary(c))
println("U:", U)
println()

include("as_mod.jl")
m, x, z = build_model(;relax_x=false, relax_z=false)
set_optimizer(m, Gurobi.Optimizer)
# set_optimizer_attributes(m, "MIPGap" => 2e-2, "TimeLimit" => 3600)
# """
# Some useful parameters for the Gurobi solver:
#     SolutionLimit = k : the search is terminated after k feasible solutions has been found
#     MIPGap = r : the search is terminated when  | best node - best integer | < r * | best node |
#     MIPGapAbs = r : the search is terminated when  | best node - best integer | < r
#     TimeLimit = t : limits the total time expended to t seconds
#     DisplayInterval = t : log lines are printed every t seconds
# See http://www.gurobi.com/documentation/8.1/refman/parameters.html for a
# complete list of valid parameters
# """
if arg == "1a"
    println("### 1a (i)")
    optimize!(m)
    obj_i = objective_value(m)
    time_i = solve_time(m)
    
    
    println("### 1a (ii)")
    unset_binary.(x)
    optimize!(m)
    obj_ii = objective_value(m)
    time_ii = solve_time(m)
    
    println("### 1a (iii)")
    # unset_binary.(x)
    unset_binary.(z)
    optimize!(m)
    obj_iii = objective_value(m)
    time_iii = solve_time(m)
    
    println("time_i = $time_i")
    println("time_ii = $time_ii")
    println("time_iii = $time_iii")
    println("obj_i = $obj_i, obj_ii = $obj_ii, obj_iii = $obj_iii")
    println("obj_i - obj_ii = ", obj_i - obj_ii)
    println("obj_ii - obj_iii = ", obj_ii - obj_iii)
elseif arg == "1b"
    println("### 1b (i)")
    optimize!(m)
    obj_i = objective_value(m)
    time_i = solve_time(m)
    
    println("### 1b (ii)")
    unset_binary.(x)
    unset_binary.(z)
    optimize!(m)
    obj_ii = objective_value(m)
    time_ii = solve_time(m)
    
    println("### 1b (iii)")
    add_cut_to_small(m)
    optimize!(m)
    obj_iii = objective_value(m)
    time_iii = solve_time(m)

    println("time_i = $time_i")
    println("time_ii = $time_ii")
    println("time_iii = $time_iii")
    println("obj_i = $obj_i, obj_ii = $obj_ii, obj_iii = $obj_iii")
    println("obj_i - obj_ii = ", obj_i - obj_ii)
    println("obj_ii - obj_iii = ", obj_ii - obj_iii)
    println("obj_i - obj_iii = ", obj_i - obj_iii)
end
"""
Some useful output & functions
"""
# obj_ip = objective_value(m)
# unset_binary.(x)
# unset_binary.(z)
# optimize!(m)
# obj_lp = objective_value(m)
# println("obj_ip = $obj_ip, obj_lp = $obj_lp, gap = $(obj_ip-obj_lp) ")

# println(solve_time(m))

x_val = sparse(value.(x.data))
z_val = sparse(value.(z))

println("x  = ")
println(x_val)
println("z = ")
println(z_val)
# println("raw x", value.(x.data))

#add_cut_to_small(m)

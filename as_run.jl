#!/usr/bin/julia
using JuMP
#using Cbc
using Gurobi
using SparseArrays
using Plots
using Printf

arg = length(ARGS) > 0 ? ARGS[1] : ""

if arg == "1a" || startswith(arg, "2")
    include("as_dat_large.jl")
else
    include("as_dat_small.jl")
end

if arg == "1a"
    global T = 125
end

println("# Data:")
println("Components:", Components)
println("T:", T)
println("d:", d)
println("c:", summary(c))
println("U:", U)
println()

if startswith(arg, "3")
    include("as_mod3.jl")
else
    include("as_mod.jl")
end
# m, x, z = build_model(;relax_x=false, relax_z=false)
# set_optimizer(m, Gurobi.Optimizer)
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

function plot_x(x, i)
    x_val = sparse(value.(x.data))
    rows, cols = findnz(x_val)
    scatter(cols, rows)
    println("x_val_$(arg)_$(i).png")
    savefig("x_val_$(arg)_$(i).png")
end
    
if startswith(arg, "1")
    global m, x, z = build_model(;relax_x=false, relax_z=false)
    set_optimizer(m, Gurobi.Optimizer)
    set_silent(m)

    if arg == "1a"
        println("### 1a (i)")
        set_silent(m)
        optimize!(m)
        obj_i = objective_value(m)
        time_i = solve_time(m)
        plot_x(x, 1)
    
        println("### 1a (ii)")
        unset_binary.(x)
        optimize!(m)
        obj_ii = objective_value(m)
        time_ii = solve_time(m)
        plot_x(x, 2)
    
        println("### 1a (iii)")
        # unset_binary.(x)
        unset_binary.(z)
        optimize!(m)
        obj_iii = objective_value(m)
        time_iii = solve_time(m)
        plot_x(x, 3)
    
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
        plot_x(x, 1)
    
        println("### 1b (ii)")
        unset_binary.(x)
        unset_binary.(z)
        optimize!(m)
        obj_ii = objective_value(m)
        time_ii = solve_time(m)
        plot_x(x, 2)
    
        println("### 1b (iii)")
        add_cut_to_small(m)
        optimize!(m)
        obj_iii = objective_value(m)
        time_iii = solve_time(m)
        plot_x(x, 3)

        println("time_i = $time_i")
        println("time_ii = $time_ii")
        println("time_iii = $time_iii")
        println("obj_i = $obj_i, obj_ii = $obj_ii, obj_iii = $obj_iii")
        println("obj_i - obj_ii = ", obj_i - obj_ii)
        println("obj_ii - obj_iii = ", obj_ii - obj_iii)
        println("obj_i - obj_iii = ", obj_i - obj_iii)
    end

elseif startswith(arg, "2") || startswith(arg, "3")
    println("### $(arg)")
    t_vals = Float64[]
    T_range = (arg == "2b" || arg == "3b") ? (50:10:700) : (50:5:200)
    
    # Sava times to file whilst doing the calculations
    open("$(arg)_time_data.txt", "w") do io
        @printf(io, "#T, time")
        for (T_idx, T_val) in enumerate(T_range)
            # Update globals. d and c must be updated based on new T
            global T = T_val
            global d = ones(1,T)*20      #cost of a maintenance occasion
            global c = [34 25 14 21 16  3 10  5  7 10]'*ones(1,T)     #costs of new components
    
            # Build model and optimize
            global m, x, z = build_model(;relax_x=false, relax_z=false)
            set_optimizer(m, Gurobi.Optimizer)
            set_optimizer_attributes(m, "MIPGap" => 2e-2, "TimeLimit" => 3600)
            set_silent(m)
            unset_binary.(x)
            if arg == "2b"
                println("Relaxing integer requirement for z")
                unset_binary.(z)
            end
            println("# T = $T:")
            optimize!(m)

            obj_i = objective_value(m)
            time_i = solve_time(m)
            push!(t_vals, time_i)
            # push!(T_vals, T)
            println("t: $time_i")

            # Print to file
            @printf(io, "%.2f, %.2f\n", T_val, time_i)
        end
    end # open(...)
    # Save plot to a .png
    plot(T_range, t_vals, xlabel="T", ylabel="Time (s)", show=true, yscale=:log10, legend=:none)
    savefig("$(arg)_time.png")
    

    # elseif arg == "2b"
    #     println("### 2b")
    #     set_silent(m)
    #     t = zeros(Float64, 14)
    #     t_x = 50:50:700
    #     plot(t_x, t .+ rand(14))
    #     # plot(rand(10))
    #     savefig("test_fig.png")
    #     optimize!(m)
    #     for i in 1:14
    #         T = i*50
    #         # Update d and c to use new T
    #         d = ones(1,T)*20      #cost of a maintenance occasion
    #         c = [34 25 14 21 16  3 10  5  7 10]'*ones(1,T)     #costs of new components
        
    #         m, x, z = build_model(;relax_x=false, relax_z=false)
    #         set_optimizer(m, Gurobi.Optimizer)
    #         set_optimizer_attributes(m, "MIPGap" => 2e-2, "TimeLimit" => 3600)
    #         set_silent(m)
    #         unset_binary.(x)
    #         unset_binary.(z)
    #         println("# T = TT:")
    #         optimize!(m)
    #         obj_i = objective_value(m)
    #         time_i = solve_time(m)
    #         t[i-0] = time_i
    #         println("t: $time_i")
    #     end
    #     plot(t_x, t, title="hello", xlabel="T", ylabel="Time (s)", show=true)
    #     savefig("2b_time.png")
    # end
else
    global m, x, z = build_model(;relax_x=false, relax_z=false)
    set_optimizer(m, Gurobi.Optimizer)
    optimize!(m)
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


# heatmap(Array(x_val))

println(x_val.colptr)

println("x  = ")
println(x_val)
println("z = ")
println(z_val)
# println("raw x", value.(x.data))

#add_cut_to_small(m)

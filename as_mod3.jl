"""
  Construct and return the model described by https://www.sciencedirect.com/science/article/pii/S0360835214000527?via%3Dihub
"""
function build_model(;relax_x::Bool = false, relax_z::Bool = false)
  # Components - the set of components
  # T - the number of time steps in the model
  # d[1,..,T] - cost of a maintenance occasion
  # c[Components, 1,..,T] - costs of new components
  # U[Components] - lives of new components
  m = Model()
  # x_ist = 1 if component i is replaced at time s and t, but not in-between, 0 otherwise.
  # z_t = 1 if maintenance occurs at time t (same as the first model)
  @variable(m, x[Components, 0:T, 1:T+1] >= 0, Bin)
  @variable(m, z[1:T] <= 1, Bin)

  # TODO: Optimize by doing: t in s+1:T+1
  # TOTHINK: Is maintenance at t=0 free?
  cost = @objective(m, Min,
    sum(c[i, s]*x[i, s, t] for i in Components, s in 1:T, t in 1:(T+1)) +
    sum(d[t]*z[t] for t in 1:T))

  CorrectIntervals = @constraint(m,
    [i in Components, t in 1:T+1],
    sum(x[i, s, t] for s in 0:t-1) == sum(x[i, t, r] for r in t+1:T+1) 
  )
    
  ReplaceWithinLife = @constraint(m,
    [i in Components, ell in 0:(T-U[i]); T >= U[i]],
    sum(x[i,0,t] for t in (ell .+ (1:U[i]))) == 1)

  ReplaceOnlyAtMaintenance = @constraint(m, [i in Components, t in 1:T],
  sum(x[i,s,t] for s in 0:t-1) <= z[t])

  return m, x, z
end

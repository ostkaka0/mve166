import OffsetArrays
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

  # We replace c with cc, such that
  #   cc[i,s,t] := T(max(d) + n*max(c)) + 1 when t-s > U[i]
  #   cc[i,s,t] := c[i,t] when t-s <= U[i] and t != T+1,
  #   0 otherwise
  # where n is last(Components) = length(Components)
  n = last(Components)
  @assert last(Components) == length(Components)
  # We have to use an OffsetArray to allow zero-indexing on 2nd argument.
  cc = OffsetArray(zeros(n, T+1, T+1), Components, 0:T, 1:T+1)
  # cc = OffsetArray(zeros(Components, T+1, T+1), 1:Components, 0:T, 1:T+1)

  for i in Components
    for s in 0:T
      for t in s+1:T+1
        if t-s > U[i]
          cc[i, s, t] = T*(maximum(d) + n*maximum(c)) + 1
        elseif t != T+1
          # Note that c[i, T+1] doesn't exist, so we let cc[i, s, T+1] it stay at zero through the the ifelse above
          cc[i, s, t] = c[i, t]
        end
      end
    end
  end
      

  # Outdated comment: # s in 1:T instead of 0:T because maintenance is free(and mandatory) at time 0. The paper had 0:T however.
  cost = @objective(m, Min,
    sum(cc[i, s, t]*x[i, s, t] for i in Components, s in 0:T, t in s+1:T+1) +
    sum(d[t]*z[t] for t in 1:T))

  ReplaceOnlyAtMaintenance = @constraint(m,
    [i in Components, t in 1:T],
    sum(x[i,s,t] for s in 0:t-1) <= z[t]
  )
  MatchingIntervals = @constraint(m,
    [i in Components, t in 1:T],
    sum(x[i, s, t] for s in 0:t-1) == sum(x[i, t, r] for r in t+1:T+1) 
  )
  MustStartMaintained = @constraint(m,
    [i in Components],
    sum(x[i,0,t] for t in 1:T+1) == 1
  )
  # # This one is not in the paper, but we must however add it because we require components to be replaced within their lives.
  # # We do this by asserting that no interval longer than U[i] can exist
  # ReplaceWithinLife = @constraint(m,
  #   [i in Components, s = 0:(T-U[i]); T >= U[i]],
  #   sum(x[i,s,t] for t in (s+U[i]+1):(T+1)) == 0
  # )
  

  return m, x, z
end

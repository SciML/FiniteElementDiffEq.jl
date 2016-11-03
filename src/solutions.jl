"""
`FEMSolution`

Holds the data for the solution to a finite element problem.

### Fields

* `fem_mesh::FEMmesh`: The finite element mesh the problem was solved on.
* `u::Array{Float64}`: The solution (at the final timepoint)
* `trueknown::Bool`: Boolean flag for if the true solution is given.
* `u_analytic::AbstractArrayOrVoid`: The true solution at the final timepoint.
* `errors`: A dictionary of the error calculations.
* `appxtrue::Bool`: Boolean flag for if u_analytic was an approximation.
* `timeseries`::AbstractArrayOrVoid`: u over time. Only saved if `save_timeseries=true`
  is specified in the solver.
* `t::AbstractArrayOrVoid`: All the t's in the solution. Only saved if `save_timeseries=true`
  is specified in the solver.
* `prob::DEProblem`: Holds the problem object used to define the problem.
* `save_timeseries::Bool`: True if solver saved the extra timepoints.

"""
type FEMSolution <: AbstractFEMSolution
  fem_mesh::FEMmesh
  u#::Array{Number}
  trueknown::Bool
  u_analytic#::AbstractArrayOrVoid
  errors#::Dict{String,Float64}
  appxtrue::Bool
  timeseries#::GrowableArray
  t#::AbstractArrayOrVoid
  prob::DEProblem
  save_timeseries::Bool
  tslocation::Int # Used in the iterator for animation
  function FEMSolution(fem_mesh::FEMmesh,u,u_analytic,sol,Du,timeSeries,t,prob;save_timeseries=true)
    errors = Dict(:L2=>getL2error(fem_mesh,sol,u),:H1=>getH1error(fem_mesh,Du,u),
                  :l∞=> maximum(abs.(u-u_analytic)), :l2=> norm(u-u_analytic,2))
    return(new(fem_mesh,u,true,u_analytic,errors,false,timeSeries,t,prob,true,0))
  end
  FEMSolution(fem_mesh,u,u_analytic,sol,Du,prob) = FEMSolution(fem_mesh::FEMmesh,u,u_analytic,sol,Du,[],[],prob,save_timeseries=false)
  function FEMSolution(fem_mesh::FEMmesh,u::AbstractArray,prob)
    return(FEMSolution(fem_mesh,u,[],[],prob,save_timeseries=false))
  end
  function FEMSolution(fem_mesh::FEMmesh,u::AbstractArray,timeseries,t,prob;save_timeseries=true)
    return(new(fem_mesh,u,false,nothing,Dict{String,Float64},false,timeseries,t,prob,save_timeseries,0))
  end
end

length(sol::FEMSolution) = length(sol.timeseries[1])

function start(sol::FEMSolution)
  sol.tslocation = 1
end

function next(sol::FEMSolution,state)
  state += 1
  sol.tslocation = state
  (sol,state)
end

function done(sol::FEMSolution,state)
  state >= length(sol)
end

eltype(sol::FEMSolution) = typeof(sol.timeseries[1][1])

"""
`S = FEMSolutionTS(timeseries::Vector{uType},numvars::Int)``
`S[i][j]` => Variable i at time j.
"""
function FEMSolutionTS{uType<:AbstractArray}(timeseries::Vector{uType},numvars::Int)
  G = Vector{typeof(timeseries[1][:,1])}(0)
  push!(G,timeseries[1][:,1])
  for j = 2:length(timeseries)
    push!(G,timeseries[j][:,1])
  end
  component_timeseries = Vector{typeof(G)}(0)
  push!(component_timeseries,G)

  if numvars > 1
    for i=2:numvars
      G = Vector{typeof(timeseries[1][:,1])}(0)
      for j = 1:length(timeseries)
        push!(G,timeseries[j][:,i])
      end
      push!(component_timeseries,G)
    end
  end
  return(component_timeseries)
end

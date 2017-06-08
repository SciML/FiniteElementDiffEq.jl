__precompile__()

module FiniteElementDiffEq

  using Reexport
  @reexport using DiffEqBase
  @reexport using DiffEqPDEBase
  import DiffEqBase: solve
  using IterativeSolvers, GenericSVD, NLsolve, Parameters,
        ChunkedArrays, Juno, VectorizedRoutines.Matlab

  include("algorithms.jl")
  include("fem_solve.jl")
  include("fem_heatintegrators.jl")


  # General Functions
  export solve

  export FEMDiffEqPoisson

  export FEMDiffEqHeatEuler, FEMDiffEqHeatImplicitEuler, FEMDiffEqHeatCrankNicholson,
         FEMDiffEqHeatSemiImplicitEuler, FEMDiffEqHeatSemiImplicitCrankNicholson

end # module

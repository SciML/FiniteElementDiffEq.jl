__precompile__()

module FiniteElementDiffEq

  using DiffEqBase, DiffEqPDEBase
  import DiffEqBase: solve
  using IterativeSolvers, GenericSVD, ForwardDiff, NLsolve, Parameters,
        ChunkedArrays, Juno, VectorizedRoutines.Matlab

  include("fem_solve.jl")
  include("fem_heatintegrators.jl")

  # General Functions
  export solve


end # module

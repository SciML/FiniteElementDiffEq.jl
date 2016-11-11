module FiniteElementDiffEq

  using DiffEqBase
  import DiffEqBase: solve
  using IterativeSolvers, Parameters, GenericSVD, ForwardDiff, NLsolve,
        ChunkedArrays, VectorizedRoutines.Matlab, RecipesBase, Juno

  import Base: size, length, start, next, done, eltype

  macro def(name, definition)
      return quote
          macro $name()
              esc($(Expr(:quote, definition)))
          end
      end
  end

  include("mesh_tools/meshes.jl")
  include("utils.jl")
  include("problems.jl")
  include("solutions.jl")
  include("plotrecipes.jl")
  include("solvers/fem_assembly.jl")
  include("mesh_tools/fem_boundary.jl")
  include("solvers/fem_error.jl")
  include("solvers/fem_solve.jl")
  include("solvers/fem_heatintegrators.jl")

  # Types
  export HeatProblem, PoissonProblem, FEMSolution, FEMmesh, SimpleMesh, FDMMesh

  # General Functions
  export solve, animate

  #FEM Functions
  export  assemblematrix, findboundary, setboundary, findbdtype, getL2error, quadpts, getH1error,
          gradu, gradbasis, quadfbasis, fem_squaremesh, CFLμ, CFLν,
          meshgrid, notime_squaremesh, parabolic_squaremesh, quadpts1
end # module

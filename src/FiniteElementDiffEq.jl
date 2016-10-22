module FiniteElementDiffEq

using DiffEqBase
import DiffEqBase: solve
using IterativeSolvers, Parameters, GenericSVD, ForwardDiff, NLsolve,
      ChunkedArrays, InplaceOps, Ranges, VectorizedRoutines.Matlab, RecipesBase

import Plots: plot
import Base: size

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
  include("premades/premade_meshes.jl")
  include("premades/premade_problems.jl")
  include("solvers/fem_solve.jl")
  include("solvers/fem_heatintegrators.jl")

  # Types
  export HeatProblem, PoissonProblem, FEMSolution, FEMmesh, SimpleMesh, FDMMesh

  # General Functions
  export solve, animate

  #FEM Example Problems
  export  prob_femheat_moving, prob_femheat_pure, prob_femheat_diffuse,
          prob_poisson_wave, prob_poisson_noisywave, prob_femheat_birthdeath,
          prob_poisson_birthdeath, prob_femheat_stochasticbirthdeath,
          prob_stokes_homogenous, prob_stokes_dirichletzero, prob_poisson_birthdeathsystem,
          prob_poisson_birthdeathinteractingsystem, prob_femheat_birthdeathinteractingsystem,
          prob_femheat_birthdeathsystem, prob_femheat_diffusionconstants,
          heatProblemExample_grayscott,heatProblemExample_gierermeinhardt

  #Example Meshes
  export  meshExample_bunny, meshExample_flowpastcylindermesh, meshExample_lakemesh,
          meshExample_Lshapemesh, meshExample_Lshapeunstructure, meshExample_oilpump,
          meshExample_wavymesh, meshExample_wavyperturbmesh

  #FEM Functions
  export  assemblematrix, findboundary, setboundary, findbdtype, getL2error, quadpts, getH1error,
          gradu, gradbasis, quadfbasis, fem_squaremesh, CFLμ, CFLν,
          meshgrid, notime_squaremesh, parabolic_squaremesh, quadpts1
end # module

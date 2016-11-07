#######
##FEM Stochastic Heat Animation Test
#######

#Generates an animation for a solution of the heat equation
#Uses Plots.jl, requires matplotlib >=1.5
using FiniteElementDiffEq, Plots#, ImageMagick
T = 5
dx = 1//2^(3)
dt = 1//2^(5)
fem_mesh = parabolic_squaremesh([0 1 0 1],dx,dt,T,:neumann)
prob = prob_femheat_stochasticbirthdeath

sol = solve(fem_mesh::FEMmesh,prob::HeatProblem,alg=:SemiImplicitCrankNicholson,save_timeseries=true,solver=:LU)

println("Generating Animation")
TEST_PLOT && animate(sol::FEMSolution;zlims=(0,3),cbar=false)

# Returns true if nothing error'd
true

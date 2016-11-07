######
##FEM Heat dx Convergence Tests
######
using FiniteElementDiffEq, DiffEqDevTools,Plots

#Travis CI Test Setting
#Not good plots, but quick for unit tests
dxs = 1.//2.^(2:-1:1)
dts = 1//2^(6) * ones(dxs) #Run at 2^-7 for best plot
#=
# Use this setup to get good plots
dt = 1//2^(14) #Small dt for Euler stability, but takes long
N = 4
topdx = 7
=#

prob = prob_femheat_moving

alg=:Euler; println(alg)
sim = test_convergence(dts::AbstractArray,dxs::AbstractArray,prob::HeatProblem,dxs;alg=alg)

alg=:ImplicitEuler; println(alg)
sim2 = test_convergence(dts::AbstractArray,dxs::AbstractArray,prob::HeatProblem,dxs;alg=alg)

alg=:CrankNicholson; println(alg)
sim3 = test_convergence(dts::AbstractArray,dxs::AbstractArray,prob::HeatProblem,dxs;alg=alg)

TEST_PLOT && plot(plot(sim),plot(sim2),plot(sim3),layout=@layout([a b c]),size=(1200,400))

#Returns true if all converge approximately dx^2
minimum([sim.𝒪est[:L2],sim2.𝒪est[:L2],sim3.𝒪est[:L2]] - 2 .<.1)

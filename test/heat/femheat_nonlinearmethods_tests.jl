######
##FEM Heat Nonlinear Test
######
using FiniteElementDiffEq, DiffEqProblemLibrary, DiffEqPDEBase

#Define a parabolic problem
prob = prob_femheat_birthdeath


#Solve it with a bunch of different algorithms, plot solution
println("Euler")
sol = solve(prob,FEMDiffEqHeatEuler())

println("Semi-implicit Euler")
sol = solve(prob,FEMDiffEqHeatSemiImplicitEuler())

println("Semi-implicit Crank Nicholson")
sol = solve(prob,FEMDiffEqHeatSemiImplicitCrankNicholson())

dx = 1//2^(2)
dt = 1//2^(3)
T = 1
fem_mesh = parabolic_squaremesh([0 1 0 1],dx,dt,T,:neumann)
println("Implicit Euler")
sol = solve(prob,FEMDiffEqHeatImplicitEuler(),autodiff=true)
println("Implicit Euler Autodiff False")
sol = solve(prob,FEMDiffEqHeatImplicitEuler(),autodiff=false)
#plot(sol)

#Returns true if nonlinear solver is correct
bool1 = maximum(abs.(sol[end] .- .777))<.01

### Stochastic Tests

#Define a parabolic problem
T = 1
dx = 1//2^(3)
dt = 1//2^(7)
fem_mesh = parabolic_squaremesh([0 1 0 1],dx,dt,T,:neumann)
prob = prob_femheat_stochasticbirthdeath


#Solve it with a bunch of different algorithms, plot solution
println("Euler")
sol = solve(prob,FEMDiffEqHeatEuler())

println("Semi-implicit Euler")
sol = solve(prob,FEMDiffEqHeatSemiImplicitEuler())

#=
# CG and GMRES require size 1 vector, breaks with numvars change
println("Semi-implicit Crank Nicholson")
sol = solve(prob,alg=:SemiImplicitCrankNicholson,solver=:CG)

println("Semi-implicit Crank Nicholson GMRES")
sol = solve(prob,alg=:SemiImplicitCrankNicholson,solver=:GMRES)
=#

#Define a quicker problem
prob =  prob_femheat_stochasticbirthdeath_fast
println("Implicit Euler")
sol = solve(prob,FEMDiffEqHeatImplicitEuler(),autodiff=true)
sol = solve(prob,FEMDiffEqHeatImplicitEuler(),autodiff=false)
#plot(sol)

bool1

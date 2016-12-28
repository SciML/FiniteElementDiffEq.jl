######
##FEM Heat Nonlinear Test
######
using FiniteElementDiffEq

#Define a parabolic problem
prob = prob_femheat_birthdeathinteractingsystem
#@code_warntype solve(prob,alg=:Euler)
sol = solve(prob,FEMDiffEqHeatEuler())

#plot(sol,plot_analytic=false,zlim=(0,2),cbar=false)

prob = prob_femheat_birthdeathsystem
sol = solve(prob,FEMDiffEqHeatImplicitEuler(),iterations=1000,autodiff=true)

#plot(sol,plot_analytic=false,zlim=(0,2),cbar=false)

true

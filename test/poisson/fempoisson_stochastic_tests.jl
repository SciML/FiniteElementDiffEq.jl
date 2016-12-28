######
##FEM Stochastic Poisson Method Tests
######
using FiniteElementDiffEq

prob = prob_poisson_noisywave

sol = solve(prob,FEMDiffEqPoisson())

#plot(sol,title=["True Deterministic Solution" "Stochastic Solution"],plot_analytic=true)
#This condition should be true with really high probability
var(sol[end]) < 8e-4

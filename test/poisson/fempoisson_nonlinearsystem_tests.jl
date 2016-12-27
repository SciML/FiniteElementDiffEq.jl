######
##FEM Poisson Nonlinear System Tests
######
using FiniteElementDiffEq


prob = prob_poisson_birthdeathsystem

sol = solve(prob)

TEST_PLOT && plot(sol,plot_analytic=false,zlim=(0,2))

#Returns true if computed solution is homogenous near 2
bool1 = maximum(abs.(sol[end] .- [2 1]))< 1e-8

### Harder system

prob = prob_poisson_birthdeathinteractingsystem

sol = solve(prob)

TEST_PLOT && plot(sol,plot_analytic=false,zlim=(0,2),cbar=false)

bool2 = maximum(abs.(sol[end] .- [2 1]))< 1e-8

bool1 && bool2

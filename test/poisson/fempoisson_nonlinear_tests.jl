######
##FEM Poisson Nonlinear Tests
######
using FiniteElementDiffEq

prob = prob_poisson_birthdeath

sol = solve(prob)

TEST_PLOT && plot(sol,plot_analytic=false,zlim=(0,2))

#Returns true if computed solution is homogenous near 2
maximum(abs.(sol.u - 2))< 1e-9

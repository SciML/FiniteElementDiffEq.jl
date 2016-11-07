######
##FEM Poisson dx Convergence Tests
######
using FiniteElementDiffEq, DiffEqDevTools#,LaTeXStrings

dxs = 1.//2.^(4:-1:2) # 4 for testing, use 7 for good graph
prob = prob_poisson_wave

sim = test_convergence(dxs::AbstractArray,prob::PoissonProblem)

#Plot Result
TEST_PLOT && plot(sim,xguide="Delta x")

#Returns true if convergence is like dx^2 in L2
sim.𝒪est[:L2]-2 <.1

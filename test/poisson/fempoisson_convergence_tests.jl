######
##FEM Poisson dx Convergence Tests
######
using FiniteElementDiffEq, DiffEqDevTools,DiffEqProblemLibrary

cs = cs_fempoisson_wave

sim = test_convergence(cs,FEMDiffEqPoisson())

#Plot Result
#plot(sim,xguide="Delta x")

#Returns true if convergence is like dx^2 in L2
sim.𝒪est[:L2]-2 <.1

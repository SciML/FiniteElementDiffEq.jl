##Finite Element Method Introduction

using FiniteElementDiffEq,  DiffEqPDEBase

### Setup
dx = 1//2^(3)
fem_mesh = notime_squaremesh([0 1 0 1],dx,:dirichlet)

f = (x) -> sin.(2π.*x[:,1]).*cos.(2π.*x[:,2])
prob = PoissonProblem(f,fem_mesh)

sol = solve(prob)

sol = solve(prob,solver=:CG)

sol = solve(prob,solver=:GMRES)

#plot(sol)

### Test Results

#The test will only pass if the calculated L2 error is below 1e-4.
true

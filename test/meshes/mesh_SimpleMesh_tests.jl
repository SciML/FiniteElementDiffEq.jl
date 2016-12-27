using FiniteElementDiffEq

pdeprob = prob_poisson_wave

res = solve(pdeprob)

mesh = SimpleMesh(fem_mesh.node,fem_mesh.elem)

true

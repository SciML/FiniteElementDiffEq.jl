##Boundary Setting Tests

using FiniteElementDiffEq, JLD

mesh = meshExample_lakemesh()

findboundary(mesh,ones(size(vec(mesh.elem))))
setboundary(mesh,:robin)

true

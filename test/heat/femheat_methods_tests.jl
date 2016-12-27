######
##FEM Heat Method Tests
######
using FiniteElementDiffEq, DiffEqDevTools

#Define a parabolic problem
prob = prob_femheat_moving7

#Solve it with a bunch of different algorithms, plot solution
println("Euler")
sol = solve(prob,alg=:Euler)
T=1 
dt = 1//2^(4) #Make faster for tests
fem_mesh = parabolic_squaremesh([0 1 0 1],dx,dt,T,:dirichlet)
println("Direct")
sol = solve(prob,alg=:ImplicitEuler,solver=:Direct)

println("LU")
sol = solve(prob,alg=:ImplicitEuler,solver=:LU)

println("QR")
sol = solve(prob,alg=:ImplicitEuler,solver=:QR)

println("SVD")
sol = solve(prob,alg=:ImplicitEuler,solver=:SVD)

println("Direct")
sol = solve(prob,alg=:CrankNicholson,solver=:Direct)

println("Cholesky")
sol = solve(prob,alg=:CrankNicholson,solver=:Cholesky)

println("CG")
sol = solve(prob,alg=:CrankNicholson,solver=:CG)

println("CG")
sol = solve(prob,alg=:ImplicitEuler,solver=:CG)

println("GMRES")
sol = solve(prob,alg=:CrankNicholson,solver=:GMRES)

#Define another parabolic problem
prob = prob_femheat_diffuse #also try heatProblemExample_pure() or heatProblemExample_diffuse()

#Solve it with a bunch of different algorithms, plot solution
println("Euler")
sol = solve(prob,alg=:Euler)

#Define a different parabolic problem
prob = prob_femheat_pure

#Solve with Euler Method
println("Euler")
sol = solve(prob,alg=:Euler)

#Choose a finer mesh, solve with Euler, and add this result to the previous as
#an approximately true solution.
prob = prob_femheat_pure11
sol2 = solve(prob,alg=:Euler)
appxtrue!(sol,sol2)
TEST_PLOT && plot(sol,plot_analytic=true,cbar=false)

sol.errors[:l2]<.005 #Returns true if res solution is near the apprxTrue res2 solution

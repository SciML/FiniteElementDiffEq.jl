######
##FEM Heat dt Convergence Tests
######
using FiniteElementDiffEq, DiffEqDevTools, Plots
#Convergences estimate has not converged in this range
#Should decrease dx/dt for better estimate
N = 2 #Number of different dt to solve at, 2 for test speed
topdt = 6 # 1//2^(topdt-1) is the max dt. Small for test speed
prob = prob_femheat_moving #also try heatProblemExample_pure() or heatProblemExample_diffuse()
dts = 1.//2.^(topdt-1:-1:N)
dxs = 1//2^(5) * ones(dts) #Run at 2^-7 for best plot


alg=:Euler; println(alg) #Unstable due to μ
sim = test_convergence(dts,dxs,prob,dts;alg=alg)

alg=:ImplicitEuler; println(alg)
sim2 = test_convergence(dts,dxs,prob,dts;alg=alg)

alg=:CrankNicholson; println(alg) #Bound by spatial discretization error at low dt, decrease dx for full convergence
dxs = 1//2^(4) * ones(dts) #Run at 2^-7 for best plot
sim3 = test_convergence(dts,dxs,prob,dts;alg=alg)

#plot(plot(sim),plot(sim2),plot(sim3),layout=@layout([a b c]),size=(1200,400))
#Note: Stabilizes in H1 due to high dx-error, reduce dx and it converges further.

#Returns true if ImplicitEuler converges like dt and
#CN convergeces like >dt^2 (approaches dt^2 as dt and dx is smaller
minimum([abs(sim2.𝒪est[:L2]-1)<.3 abs(sim3.𝒪est[:L2]-2)<.1])

######
##FEM Heat dt Convergence Tests
######
using FiniteElementDiffEq, DiffEqDevTools, DiffEqProblemLibrary
#Convergences estimate has not converged in this range
#Should decrease dx/dt for better estimate

cs = cs_femheat_moving_dt
alg=:Euler; println(alg) #Unstable due to μ
sim = test_convergence(cs;alg=alg)

alg=:ImplicitEuler; println(alg)
sim2 = test_convergence(cs;alg=alg)

alg=:CrankNicholson; println(alg) #Bound by spatial discretization error at low dt, decrease dx for full convergence
cs = cs_femheat_moving_faster_dt
sim3 = test_convergence(cs;alg=alg)

#plot(plot(sim),plot(sim2),plot(sim3),layout=@layout([a b c]),size=(1200,400))
#Note: Stabilizes in H1 due to high dx-error, reduce dx and it converges further.

#Returns true if ImplicitEuler converges like dt and
#CN convergeces like >dt^2 (approaches dt^2 as dt and dx is smaller
minimum([abs(sim2.𝒪est[:L2]-1)<.3 abs(sim3.𝒪est[:L2]-2)<.1])

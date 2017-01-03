######
##FEM Heat dx Convergence Tests
######
using FiniteElementDiffEq, DiffEqDevTools

#Travis CI Test Setting
#Not good plots, but quick for unit tests
cs = cs_femheat_moving_dx

alg=FEMDiffEqHeatEuler(); println(alg) #Unstable due to μ
sim = test_convergence(cs,alg)

alg=FEMDiffEqHeatImplicitEuler(); println(alg)
sim2 = test_convergence(cs,alg)

alg=FEMDiffEqHeatCrankNicholson(); println(alg) #Bound by spatial discretization error at low dt, decrease dx for full convergence
cs = cs_femheat_moving_faster_dt
sim3 = test_convergence(cs,alg)

#plot(plot(sim),plot(sim2),plot(sim3),layout=@layout([a b c]),size=(1200,400))

#Returns true if all converge approximately dx^2
minimum([sim.𝒪est[:L2],sim2.𝒪est[:L2],sim3.𝒪est[:L2]] - 2 .<.1)

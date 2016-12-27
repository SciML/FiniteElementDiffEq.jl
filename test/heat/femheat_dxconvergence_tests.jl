######
##FEM Heat dx Convergence Tests
######
using FiniteElementDiffEq, DiffEqDevTools,Plots

#Travis CI Test Setting
#Not good plots, but quick for unit tests
cs = cs_femheat_moving_dx

alg=:Euler; println(alg)
sim = test_convergence(cs;alg=alg)

alg=:ImplicitEuler; println(alg)
sim2 = test_convergence(cs;alg=alg)

alg=:CrankNicholson; println(alg)
sim3 = test_convergence(cs;alg=alg)

TEST_PLOT && plot(plot(sim),plot(sim2),plot(sim3),layout=@layout([a b c]),size=(1200,400))

#Returns true if all converge approximately dx^2
minimum([sim.𝒪est[:L2],sim2.𝒪est[:L2],sim3.𝒪est[:L2]] - 2 .<.1)

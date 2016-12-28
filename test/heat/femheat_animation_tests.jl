#######
##FEM Heat Animation Test
#######

#Generates an animation for a solution of the heat equation
#Uses Plots.jl, requires matplotlib >=1.5
using FiniteElementDiffEq, DiffEqProblemLibrary
#using Plots, ImageMagick; gr()
prob = prob_femheat_moving

sol = solve(prob,FEMDiffEqHeatEuler(),save_timeseries=true)
println("Generating Animation")
sol.tslocation = 1
#plot(sol)
#animate(sol,"test_animation.gif";zlims=(0,.1),cbar=false)

## Should have moved off the frame.
maximum(sol[end]) .< 1e-6

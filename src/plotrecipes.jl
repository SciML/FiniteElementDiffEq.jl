#=
"""
`animate(sol::FEMSolution)`

Plots an animation of the solution. Requires `save_timeseries=true` was enabled in the solver.
"""
function animate(sol::FEMSolution;filename="tmp.gif",fps=15,kw...)
  atomloaded = isdefined(Main,:Atom)
  anim = Plots.Animation()
  for j=1:length(sol.timeseries[1])
    plot(sol,tslocation=j;kw...)
    Plots.frame(anim)
    atomloaded ? Main.Atom.progress(j/length(sol.timeseries[1])) : nothing #Use Atom's progressbar if loaded
  end
  gif(anim,filename,fps=fps)
end
=#

@recipe function f(sol::FEMSolution;plot_analytic=false)
  if sol.tslocation==0 #Plot solution at end
    out = Any[]
    for i = 1:size(sol.u,2)
      push!(out,sol.u[:,i])
    end
    if plot_analytic
      for i = 1:size(sol.u,2)
        push!(out,sol.u_analytic[:,i])
      end
    end
  else #use timeseries
    out = Any[]
    for i = 1:sol.prob.numvars
      push!(out,sol.timeseries[i][sol.tslocation])
    end
  end
  seriestype --> :surface
  layout --> length(out)
  sol.fem_mesh.node[:,1], sol.fem_mesh.node[:,2], out
end


@recipe function f(mesh::Mesh)
  seriestype --> :surface #:plot_trisurf
  #triangles  --> mesh.elem-1
  mesh.node[:,1], mesh.node[:,2], ones(mesh.node[:,1])
end

# mesh = meshExample_lakemesh()
# PyPlot.plot_trisurf(mesh.node[:,1],mesh.node[:,2],ones(mesh.node[:,2]),triangles=mesh.elem-1)

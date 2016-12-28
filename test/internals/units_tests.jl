using Unitful, FiniteElementDiffEq

### Setup
dx = (1//2^(3))u"m"
mesh = notime_squaremesh([0u"m" 1u"m" 0u"m" 1u"m"],dx,:dirichlet) #Fails at meshgrid because unit ranges not supported

f = (x) -> sin(2π.*map((y)->y.val,x[:,1])).*cos(2π.*map((y)->y.val,x[:,2]))
prob = PoissonProblem(f,mesh,D=1u"m^2/s")

#=
sol = solve(mesh,prob)
=#

#Define a parabolic problem
T = 1u"s"
dt = (1//2^(7))u"s"
mesh = parabolic_squaremesh([0u"m" 1u"m" 0u"m" 1u"m"],dx,dt,T,:dirichlet)
f = (u,x,t)  -> (float(ones(size(x,1))))u"N" - .5u
u0_func = (x) -> map((x)->(x)u"N",zeros(size(x,1)))
u0 = u0_func(mesh.node)
prob = HeatProblem(u0,f,mesh)
#=
println("Euler")
sol = solve(prob,alg=:Euler)
=#

true

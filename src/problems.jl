type HeatProblem{MeshType,F,F2,F3,F4,F5,F6,F7,DiffType} <: AbstractHeatProblem
  u0::F5
  Du::F2
  f::F
  gD::F3
  gN::F4
  analytic::F7
  knownanalytic::Bool
  islinear::Bool
  numvars::Int
  σ::F6
  stochastic::Bool
  noisetype::Symbol
  D::DiffType
  mesh::MeshType
end

function HeatProblem(analytic,Du,f,mesh;gN=nothing,σ=nothing,noisetype=:White,numvars=nothing,D=nothing)
  islinear = numparameters(f)==2
  knownanalytic = true
  u0(x) = analytic(0,x)
  numvars = size(u0([0 0
                     0 0
                     0 0]),2)
  gD = analytic
  if gN == nothing
    gN=(t,x)->zeros(size(x,1),numvars)
  end
  if σ==nothing
    stochastic=false
    σ=(t,x)->zeros(size(x,1),numvars)
  else
    stochastic=true
  end
  if D == nothing
    if numvars == 1
      D = 1.0
    else
      D = ones(1,numvars)
    end
  end

  if numvars==0
    u = copy(u0(mesh.node))
    numvars = size(u,2)
    if gD == nothing
      gD=(t,x)->zeros(size(x,1),numvars)
    end
    if gN == nothing
      gN=(t,x)->zeros(size(x,1),numvars)
    end
    if D == nothing
      if numvars == 1
        D = 1.0
      else
        D = ones(1,numvars)
      end
    end
  end

  HeatProblem(u0,Du,f,gD,gN,analytic,knownanalytic,islinear,numvars,σ,stochastic,noisetype,D,mesh)
end


function HeatProblem(u0,f,mesh;gD=nothing,gN=nothing,σ=nothing,noisetype=:White,numvars=nothing,D=nothing)
  if σ==nothing
    stochastic=false
    σ=(t,x)->zeros(size(x,1))
  else
    stochastic=true
  end
  islinear = numparameters(f)==2
  knownanalytic = false
  if islinear
    if u0==nothing
      u0=(x)->zeros(size(x,1))
    end
    if gD == nothing
      gD=(t,x)->zeros(size(x,1))
    end
    if gN == nothing
      gN=(t,x)->zeros(size(x,1))
    end
    if D == nothing
      D = 1.0
    end
    numvars = 1
  end
  if !islinear #nonlinear
    if u0==nothing && numvars == nothing
      warn("u0 and numvars must be given. numvars assumed 1.")
      numvars = 1
      u0=(x)->zeros(size(x,1),numvars)
      if gD == nothing
        gD=(t,x)->zeros(size(x,1),numvars)
      end
      if gN == nothing
        gN=(t,x)->zeros(size(x,1),numvars)
      end
      if D == nothing
        D = 1.0
      end
    elseif u0==nothing #numvars!=nothing
      u0=(x)->zeros(size(x,1),numvars) #Default to zero
      if gD == nothing
        gD=(t,x)->zeros(size(x,1),numvars)
      end
      if gN == nothing
        gN=(t,x)->zeros(size(x,1),numvars)
      end
      if D == nothing
        D = ones(1,numvars)
      end
    elseif numvars==nothing #If u0 is given but numvars is not, we're still okay. Generate from size in function.
      numvars=0 #Placeholder, update gD and gN in solver
    end
  end

  if numvars==0
    u = copy(u0(mesh.node))
    numvars = size(u,2)
    if gD == nothing
      gD=(t,x)->zeros(size(x,1),numvars)
    end
    if gN == nothing
      gN=(t,x)->zeros(size(x,1),numvars)
    end
    if D == nothing
      if numvars == 1
        D = 1.0
      else
        D = ones(1,numvars)
      end
    end
  end

  HeatProblem(u0,(x)->0,f,gD,gN,(x)->0,knownanalytic,islinear,numvars,σ,stochastic,noisetype,D,mesh)
end

type PoissonProblem{MeshType} <: AbstractPoissonProblem
  "f: Forcing function in the Poisson problem"
  f#::Function
  "analytic: Solution to the Poisson problem"
  analytic::Function
  "Du: Gradient of the solution to the Poisson problem"
  Du::Function
  "gD: dirichlet Boundary Data"
  gD#::Nullable{Function}
  "gN: neumann Boundary Data"
  gN#::Nullable{Function}
  "knownanalytic: Boolean which states whether the solution function is given"
  knownanalytic::Bool
  "islinear: Boolean which states whether the problem is linear or nonlinear"
  islinear::Bool
  u0::Function
  numvars::Int
  σ::Function
  stochastic::Bool
  noisetype::Symbol
  D#::AbstractArray
  mesh::MeshType
end

function PoissonProblem(f,analytic,Du,mesh;gN=nothing,σ=nothing,u0=nothing,noisetype=:White,numvars=nothing,D=nothing)
  gD = analytic
  numvars = size(analytic([0 0
                      0 0
                      0 0]),2)
  islinear = numparameters(f)==1
  if gN == nothing
    gN=(x)->zeros(size(x,1),numvars)
  end
  if u0==nothing
    u0=(x)->zeros(size(x,1),numvars)
  end
  if D == nothing
    if numvars == 1
      D = 1.0
    else
      D = ones(1,numvars)
    end
  end
  if σ==nothing
    stochastic=false
    σ=(x)->zeros(size(x,1),numvars)
  else
    stochastic=true
  end

  u = u0(mesh.node)

  if numvars==0
    numvars = size(u,2)
    if gD == nothing
      gD=(x)->zeros(size(x,1),numvars)
    end
    if gN == nothing
      gN=(x)->zeros(size(x,1),numvars)
    end
    if D == nothing
      if numvars == 1
        D = 1.0
      else
        D = ones(1,numvars)
      end
    end
  end
  PoissonProblem(f,analytic,Du,analytic,gN,true,islinear,u0,numvars,σ,stochastic,noisetype,D,mesh)
end
function PoissonProblem(f,mesh;gD=nothing,gN=nothing,u0=nothing,σ=nothing,noisetype=:White,numvars=nothing,D=nothing)
  if σ==nothing
    stochastic=false
    σ=(x)->zeros(size(x,1))
  else
    stochastic = true
  end
  islinear = numparameters(f)==1
  if islinear && u0==nothing
    u0=(x)->zeros(size(x,1))
    if gD == nothing
      gD=(x)->zeros(size(x,1))
    end
    if gN == nothing
      gN=(x)->zeros(size(x,1))
    end
    if D == nothing
      D = 1.0
    end
    numvars = 1
  end
  if !islinear #nonlinear
    if u0==nothing && numvars == nothing
      warn("u0 and numvars must be given. numvars assumed 1.")
      numvars = 1
      u0=(x)->zeros(size(x,1))
      if gD == nothing
        gD=(x)->zeros(size(x,1))
      end
      if gN == nothing
        gN=(x)->zeros(size(x,1))
      end
      if D == nothing
        D = 1.0
      end
    elseif u0==nothing #numvars!=nothing
      u0=(x)->zeros(size(x,1),numvars) #Default to zero
      if gD == nothing
        gD=(x)->zeros(size(x,1),numvars)
      end
      if gN == nothing
        gN=(x)->zeros(size(x,1),numvars)
      end
      if D == nothing
        D = ones(1,numvars)
      end
    elseif numvars==nothing #If u0 is given but numvars is not, we're still okay. Generate from size in function.
      numvars=0 #Placeholder, update gD and gN in solver
    end
  end

  u = u0(mesh.node)

  if numvars==0
    numvars = size(u,2)
    if gD == nothing
      gD=(x)->zeros(size(x,1),numvars)
    end
    if gN == nothing
      gN=(x)->zeros(size(x,1),numvars)
    end
    if D == nothing
      if numvars == 1
        D = 1.0
      else
        D = ones(1,numvars)
      end
    end
  end
  PoissonProblem(f,(x)->0,(x)->0,gD,gN,false,islinear,u0,numvars,σ,stochastic,noisetype,D,mesh)
end

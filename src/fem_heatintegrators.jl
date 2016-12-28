immutable FEMHeatIntegrator{T1,T2,T3,F1,F2,F3,F4,uElType,uType,nodeType,AType,tType,DiffType,MinvType,StiffType}
  N::Int
  NT::Int
  dt::tType
  t::tType
  Minv::MinvType
  D::DiffType
  A::Array{uElType,2}
  freenode::Vector{Int}
  f::F1
  gD::F2
  gN::F3
  u::uType
  node::nodeType
  elem::Array{Int,2}
  area::AType
  bdnode::Vector{Int}
  mid::Array{uElType,3}
  dirichlet::Array{Int,2}
  neumann::Array{Int,2}
  islinear::Bool
  numvars::Int
  sqrtdt::Float64
  σ::F4
  noisetype::Symbol
  numiters::Int
  save_timeseries::Bool
  timeseries::Vector{uType}
  ts::Vector{tType}
  solver::Symbol
  autodiff::Bool
  method::Symbol
  show_trace::Bool
  iterations::Int
  timeseries_steps::Int
  progressbar::Bool
  progress_steps::Int
  progressbar_name::String
end

@def femheat_footer begin
  u[bdnode] = gD(i*dt,node)[bdnode]
  if save_timeseries && i%timeseries_steps==0
    push!(timeseries,copy(u))
    push!(ts,t)
  end
  if progressbar && i%progress_steps==0
    Juno.msg(prog,"dt="*string(dt))
    Juno.progress(prog,i/numiters)
  end
end

@def femheat_deterministicimplicitlinearsolve begin
  if solver==:Direct || solver==:Cholesky || solver==:QR || solver==:LU || solver==:SVD
    u[freenode,:] = lhs\(Dinv.*rhs(i,u))
  elseif solver==:CG
    for j=1:size(u,2)
      u[freenode,j],ch = cg!(u[freenode,j],lhs,Dinv.*rhs(i,u)[:,j]) # Requires Vector, need to change rhs
    end
  elseif solver==:GMRES
    for j=1:size(u,2)
      u[freenode,j],ch = gmres!(u[freenode,j],lhs,Dinv.*rhs(i,u)[:,j]) # Requires Vector, need to change rhs
    end
  end
end

@def femheat_stochasticimplicitlinearsolve begin
  dW = next(rands)
  if solver==:Direct || solver==:Cholesky || solver==:QR || solver==:LU || solver==:SVD
    u[freenode,:] = lhs\(Dinv.*rhs(i,u,dW))
  elseif solver==:CG
    for j=1:size(u,2)
      u[freenode,j],ch = cg!(u[freenode,j],lhs,Dinv.*rhs(i,u,dW)[:,j]) # Requires Vector, need to change rhs
    end
  elseif solver==:GMRES
    for j=1:size(u,2)
      u[freenode,j],ch = gmres!(u[freenode,j],lhs,Dinv.*rhs(i,u,dW)[:,j]) # Requires Vector, need to change rhs
    end
  end
end

@def femheat_implicitpreamble begin
  Dinv = D.^(-1)
  if solver==:Cholesky
    lhs = cholfact(lhs) # Requires positive definite, may be violated
  elseif solver==:LU
    lhs = lufact(lhs)
  elseif solver==:QR
    lhs = qrfact(lhs) #More stable, slower than LU
  elseif solver==:SVD
    lhs = svdfact(lhs)
  end
end

@def femheat_nonlinearsolvestochasticloop begin
  u = vec(u)
  uOld = copy(u)
  dW = next(rands)
  nlres = NLsolve.nlsolve((u,resid)->rhs!(i,u,resid,dW,uOld),uOld,autodiff=autodiff,method=method,show_trace=show_trace,iterations=iterations)
  u = nlres.zero
  if numvars > 1
    u = reshape(u,N,numvars)
  end
end

@def femheat_nonlinearsolvedeterministicloop begin
  u = vec(u)
  uOld = copy(u)
  nlres = NLsolve.nlsolve((u,resid)->rhs!(i,u,resid,uOld),uOld,autodiff=autodiff,method=method,show_trace=show_trace,iterations=iterations)
  u = nlres.zero
  if numvars > 1
    u = reshape(u,N,numvars)
  end
end

@def femheat_nonlinearsolvepreamble begin
  uOld = similar(vec(u))
end

@def femheat_deterministicpreamble begin
  @unpack N,NT,dt,t,Minv,D,A,freenode,f,gD,gN,u,node,elem,area,bdnode,mid,dirichlet,neumann,islinear,numvars,numiters,save_timeseries,timeseries,ts,solver,autodiff,method,show_trace,iterations,timeseries_steps,progressbar,progress_steps, progressbar_name = integrator
  progressbar && (prog = Juno.ProgressBar(name=progressbar_name))
end

@def femheat_stochasticpreamble begin
  @unpack sqrtdt,σ,noisetype = integrator
  rands = getNoise(u,node,elem,noisetype=noisetype)
end

@def femheat_postamble begin
  progressbar && Juno.done(prog)
  if ts[end] != t
    push!(timeseries,copy(u))
    push!(ts,t)
  end
  u,timeseries,ts
end

function femheat_solve{F1,F2,F3,F4,uElType,uType,nodeType,AType,tType,DiffType,MinvType,StiffType}(integrator::FEMHeatIntegrator{true,FEMDiffEqHeatEuler,false,F1,F2,F3,F4,uElType,uType,nodeType,AType,tType,DiffType,MinvType,StiffType})
  @femheat_deterministicpreamble
  K = eye(N) - dt*Minv*D*A #D okay since numVar = 1 for linear
  @inbounds for i=1:numiters
    u[freenode,:] = K[freenode,freenode]*u[freenode,:] + (Minv*dt*quadfbasis((x)->f(t,x),(x)->gD(t,x),(x)->gN(t,x),
                A,u,node,elem,area,bdnode,mid,N,NT,dirichlet,neumann,islinear,numvars))[freenode,:]
    t += dt
    @femheat_footer
  end
  @femheat_postamble
end

function femheat_solve{F1,F2,F3,F4,uElType,uType,nodeType,AType,tType,DiffType,MinvType,StiffType}(integrator::FEMHeatIntegrator{true,FEMDiffEqHeatEuler,true,F1,F2,F3,F4,uElType,uType,nodeType,AType,tType,DiffType,MinvType,StiffType})
  @femheat_deterministicpreamble
  @femheat_stochasticpreamble
  K = eye(N) - dt*Minv*D*A #D okay since numVar = 1 for linear
  @inbounds for i=1:numiters
    dW = next(rands)
    u[freenode,:] = K[freenode,freenode]*u[freenode,:] + (Minv*dt*quadfbasis((x)->f(t,x),(x)->gD(t,x),(x)->gN(t,x),
                A,u,node,elem,area,bdnode,mid,N,NT,dirichlet,neumann,islinear,numvars))[freenode,:] +
                (sqrtdt.*dW.*Minv*quadfbasis((x)->σ(t,x),(x)->gD(t,x),(x)->gN(t,x),
                            A,u,node,elem,area,bdnode,mid,N,NT,dirichlet,neumann,islinear,numvars))[freenode,:]
    t += dt
    @femheat_footer
  end
  @femheat_postamble
end

function femheat_solve{F1,F2,F3,F4,uElType,uType,nodeType,AType,tType,DiffType,MinvType,StiffType}(integrator::FEMHeatIntegrator{false,FEMDiffEqHeatEuler,false,F1,F2,F3,F4,uElType,uType,nodeType,AType,tType,DiffType,MinvType,StiffType})
  @femheat_deterministicpreamble
  @inbounds for i=1:numiters
    u[freenode,:] = u[freenode,:] - D.*(dt*Minv[freenode,freenode]*A[freenode,freenode]*u[freenode,:]) + (Minv*dt*quadfbasis((x,u)->f(t,x,u),(x)->gD(t,x),(x)->gN(t,x),
            A,u,node,elem,area,bdnode,mid,N,NT,dirichlet,neumann,islinear,numvars))[freenode,:]
    t += dt
    @femheat_footer
  end
  @femheat_postamble
end

function femheat_solve{F1,F2,F3,F4,uElType,uType,nodeType,AType,tType,DiffType,MinvType,StiffType}(integrator::FEMHeatIntegrator{false,FEMDiffEqHeatEuler,true,F1,F2,F3,F4,uElType,uType,nodeType,AType,tType,DiffType,MinvType,StiffType})
  @femheat_deterministicpreamble
  @femheat_stochasticpreamble
  @inbounds for i=1:numiters
    dW = next(rands)
    u[freenode,:] = u[freenode,:] - D.*(dt*Minv[freenode,freenode]*A[freenode,freenode]*u[freenode,:]) + (Minv*dt*quadfbasis((x,u)->f(t,x,u),(x)->gD(t,x),(x)->gN(t,x),
                A,u,node,elem,area,bdnode,mid,N,NT,dirichlet,neumann,islinear,numvars))[freenode,:] +
                (sqrtdt.*dW.*(Minv*quadfbasis((x,u)->σ(t,x,u),(x)->gD(t,x),(x)->gN(t,x),
                            A,u,node,elem,area,bdnode,mid,N,NT,dirichlet,neumann,islinear,numvars)))[freenode,:]
    t += dt
    @femheat_footer
  end
  @femheat_postamble
end

function femheat_solve{F1,F2,F3,F4,uElType,uType,nodeType,AType,tType,DiffType,MinvType,StiffType}(integrator::FEMHeatIntegrator{true,FEMDiffEqHeatImplicitEuler,true,F1,F2,F3,F4,uElType,uType,nodeType,AType,tType,DiffType,MinvType,StiffType})
  @femheat_deterministicpreamble
  @femheat_stochasticpreamble
  K = eye(N) + dt*Minv*D*A #D okay since numVar = 1 for linear
  lhs = K[freenode,freenode]
  rhs(i,u,dW) = u[freenode,:] + (Minv*dt*quadfbasis((x)->f((i)*dt,x),(x)->gD((i)*dt,x),(x)->gN((i)*dt,x),A,u,node,elem,area,bdnode,mid,N,NT,dirichlet,neumann,islinear,numvars))[freenode,:] +
              (sqrtdt.*dW.*(Minv*quadfbasis((x)->σ((i-1)*dt,x),(x)->gD((i-1)*dt,x),(x)->gN((i-1)*dt,x),
                          A,u,node,elem,area,bdnode,mid,N,NT,dirichlet,neumann,islinear,numvars)))[freenode,:]
  @femheat_implicitpreamble
  @inbounds for i=1:numiters
    dW = next(rands)
    t += dt
    @femheat_stochasticimplicitlinearsolve
    @femheat_footer
  end
  @femheat_postamble
end

function femheat_solve{F1,F2,F3,F4,uElType,uType,nodeType,AType,tType,DiffType,MinvType,StiffType}(integrator::FEMHeatIntegrator{true,FEMDiffEqHeatImplicitEuler,false,F1,F2,F3,F4,uElType,uType,nodeType,AType,tType,DiffType,MinvType,StiffType})
  @femheat_deterministicpreamble
  K = eye(N) + dt*Minv*D*A #D okay since numVar = 1 for linear
  lhs = K[freenode,freenode]
  rhs(i,u) = u[freenode,:] + (Minv*dt*quadfbasis((x)->f((i)*dt,x),(x)->gD((i)*dt,x),(x)->gN((i)*dt,x),A,u,node,elem,area,bdnode,mid,N,NT,dirichlet,neumann,islinear,numvars))[freenode,:]
  @femheat_implicitpreamble
  @inbounds for i=1:numiters
    t += dt
    @femheat_deterministicimplicitlinearsolve
    @femheat_footer
  end
  @femheat_postamble
end

function femheat_solve{F1,F2,F3,F4,uElType,uType,nodeType,AType,tType,DiffType,MinvType,StiffType}(integrator::FEMHeatIntegrator{true,FEMDiffEqHeatCrankNicholson,true,F1,F2,F3,F4,uElType,uType,nodeType,AType,tType,DiffType,MinvType,StiffType})
  @femheat_deterministicpreamble
  @femheat_stochasticpreamble
  Km = eye(N) - dt*Minv*D*A/2 #D okay since numVar = 1 for linear
  Kp = eye(N) + dt*Minv*D*A/2 #D okay since numVar = 1 for linear
  lhs = Kp[freenode,freenode]
  rhs(i,u,dW) = Km[freenode,freenode]*u[freenode,:] + (Minv*dt*quadfbasis((x)->f((i-.5)*dt,x),(x)->gD((i-.5)*dt,x),(x)->gN((i-.5)*dt,x),A,u,node,elem,area,bdnode,mid,N,NT,dirichlet,neumann,islinear,numvars))[freenode,:] +
              (sqrtdt.*dW.*(Minv*quadfbasis((x)->σ((i-1)*dt,x),(x)->gD((i-1)*dt,x),(x)->gN((i-1)*dt,x),
                          A,u,node,elem,area,bdnode,mid,N,NT,dirichlet,neumann,islinear,numvars)))[freenode,:]
  @femheat_implicitpreamble
  @inbounds for i=1:numiters
    dW = next(rands)
    t += dt
    @femheat_stochasticimplicitlinearsolve
    @femheat_footer
  end
  @femheat_postamble
end

function femheat_solve{F1,F2,F3,F4,uElType,uType,nodeType,AType,tType,DiffType,MinvType,StiffType}(integrator::FEMHeatIntegrator{true,FEMDiffEqHeatCrankNicholson,false,F1,F2,F3,F4,uElType,uType,nodeType,AType,tType,DiffType,MinvType,StiffType})
  @femheat_deterministicpreamble
  Km = eye(N) - dt*Minv*D*A/2 #D okay since numVar = 1 for linear
  Kp = eye(N) + dt*Minv*D*A/2 #D okay since numVar = 1 for linear
  lhs = Kp[freenode,freenode]
  rhs(i,u) = Km[freenode,freenode]*u[freenode,:] + (Minv*dt*quadfbasis((x)->f((i-.5)*dt,x),(x)->gD((i-.5)*dt,x),(x)->gN((i-.5)*dt,x),A,u,node,elem,area,bdnode,mid,N,NT,dirichlet,neumann,islinear,numvars))[freenode,:]
  @femheat_implicitpreamble
  @inbounds for i=1:numiters
    t += dt
    @femheat_deterministicimplicitlinearsolve
    @femheat_footer
  end
  @femheat_postamble
end

function femheat_solve{F1,F2,F3,F4,uElType,uType,nodeType,AType,tType,DiffType,MinvType,StiffType}(integrator::FEMHeatIntegrator{false,FEMDiffEqHeatSemiImplicitEuler,false,F1,F2,F3,F4,uElType,uType,nodeType,AType,tType,DiffType,MinvType,StiffType}) #Incorrect for system with different diffusions
  @femheat_deterministicpreamble
  Dinv = D.^(-1)
  K = eye(N) + dt*Minv*A
  lhs = K[freenode,freenode]
  rhs(i,u) = u[freenode,:] + (Minv*dt*quadfbasis((x,u)->f((i)*dt,x,u),(x)->gD((i)*dt,x),(x)->gN((i)*dt,x),
              A,u,node,elem,area,bdnode,mid,N,NT,dirichlet,neumann,islinear,numvars))[freenode,:]
  @femheat_implicitpreamble
  @inbounds for i=1:numiters
    t += dt
    @femheat_deterministicimplicitlinearsolve
    @femheat_footer
  end
  @femheat_postamble
end

function femheat_solve{F1,F2,F3,F4,uElType,uType,nodeType,AType,tType,DiffType,MinvType,StiffType}(integrator::FEMHeatIntegrator{false,FEMDiffEqHeatSemiImplicitEuler,true,F1,F2,F3,F4,uElType,uType,nodeType,AType,tType,DiffType,MinvType,StiffType}) #Incorrect for system with different diffusions
  @femheat_deterministicpreamble
  @femheat_stochasticpreamble
  Dinv = D.^(-1)
  K = eye(N) + dt*Minv*A
  lhs = K[freenode,freenode]
  rhs(i,u,dW) = u[freenode,:] + (Minv*dt*quadfbasis((x,u)->f((i)*dt,x,u),(x)->gD((i)*dt,x),(x)->gN((i)*dt,x),
              A,u,node,elem,area,bdnode,mid,N,NT,dirichlet,neumann,islinear,numvars))[freenode,:] +
              (sqrtdt.*dW.*(Minv*quadfbasis((x,u)->σ((i-1)*dt,x,u),(x)->gD((i-1)*dt,x),(x)->gN((i-1)*dt,x),
                          A,u,node,elem,area,bdnode,mid,N,NT,dirichlet,neumann,islinear,numvars)))[freenode,:]
  @femheat_implicitpreamble
  @inbounds for i=1:numiters
    dW = next(rands)
    t += dt
    @femheat_stochasticimplicitlinearsolve
    @femheat_footer
  end
  @femheat_postamble
end

function femheat_solve{F1,F2,F3,F4,uElType,uType,nodeType,AType,tType,DiffType,MinvType,StiffType}(integrator::FEMHeatIntegrator{false,FEMDiffEqHeatSemiImplicitCrankNicholson,false,F1,F2,F3,F4,uElType,uType,nodeType,AType,tType,DiffType,MinvType,StiffType}) #Incorrect for system with different diffusions
  @femheat_deterministicpreamble
  Dinv = D.^(-1)
  Km = eye(N) - dt*Minv*A/2
  Kp = eye(N) + dt*Minv*A/2
  lhs = Kp[freenode,freenode]
  rhs(i,u) = Km[freenode,freenode]*u[freenode,:] + (Minv*dt*quadfbasis((x,u)->f((i-.5)*dt,x,u),(x)->gD((i-.5)*dt,x),(x)->gN((i-.5)*dt,x),
              A,u,node,elem,area,bdnode,mid,N,NT,dirichlet,neumann,islinear,numvars))[freenode,:]
  @femheat_implicitpreamble
  @inbounds for i=1:numiters
    t += dt
    @femheat_deterministicimplicitlinearsolve
    @femheat_footer
  end
  @femheat_postamble
end

function femheat_solve{F1,F2,F3,F4,uElType,uType,nodeType,AType,tType,DiffType,MinvType,StiffType}(integrator::FEMHeatIntegrator{false,FEMDiffEqHeatSemiImplicitCrankNicholson,true,F1,F2,F3,F4,uElType,uType,nodeType,AType,tType,DiffType,MinvType,StiffType}) #Incorrect for system with different diffusions
  @femheat_deterministicpreamble
  @femheat_stochasticpreamble
  Dinv = D.^(-1)
  Km = eye(N) - dt*Minv*A/2
  Kp = eye(N) + dt*Minv*A/2
  lhs = Kp[freenode,freenode]
  rhs(i,u,dW) = Km[freenode,freenode]*u[freenode,:] + (Minv*dt*quadfbasis((x,u)->f((i-.5)*dt,x,u),(x)->gD((i-.5)*dt,x),(x)->gN((i-.5)*dt,x),
              A,u,node,elem,area,bdnode,mid,N,NT,dirichlet,neumann,islinear,numvars))[freenode,:] +
              (sqrtdt.*dW.*(Minv*quadfbasis((x,u)->σ((i-1)*dt,x,u),(x)->gD((i-1)*dt,x),(x)->gN((i-1)*dt,x),
                          A,u,node,elem,area,bdnode,mid,N,NT,dirichlet,neumann,islinear,numvars)))[freenode,:]
  @femheat_implicitpreamble
  @inbounds for i=1:numiters
    dW = next(rands)
    t += dt
    @femheat_stochasticimplicitlinearsolve
    @femheat_footer
  end
  @femheat_postamble
end

function femheat_solve{F1,F2,F3,F4,uElType,uType,nodeType,AType,tType,DiffType,MinvType,StiffType}(integrator::FEMHeatIntegrator{false,FEMDiffEqHeatImplicitEuler,false,F1,F2,F3,F4,uElType,uType,nodeType,AType,tType,DiffType,MinvType,StiffType})
  @femheat_deterministicpreamble
  function rhs!(i,u,resid,uOld)
    u = reshape(u,N,numvars)
    uOld = reshape(uOld,N,numvars)
    resid = reshape(resid,N,numvars)
    resid[freenode,:] = u[freenode,:] - uOld[freenode,:] + D.*(dt*Minv[freenode,freenode]*A[freenode,freenode]*u[freenode,:]) -
    (Minv*dt*quadfbasis((x,u)->f((i)*dt,x,u),(x)->gD((i)*dt,x),(x)->gN((i)*dt,x),A,u,node,elem,area,bdnode,mid,N,NT,dirichlet,neumann,islinear,numvars))[freenode,:]
    u = vec(u)
    resid = vec(resid)
  end
  @femheat_nonlinearsolvepreamble
  @inbounds for i=1:numiters
    t += dt
    @femheat_nonlinearsolvedeterministicloop
    @femheat_footer
  end
  @femheat_postamble
end

function femheat_solve{F1,F2,F3,F4,uElType,uType,nodeType,AType,tType,DiffType,MinvType,StiffType}(integrator::FEMHeatIntegrator{false,FEMDiffEqHeatImplicitEuler,true,F1,F2,F3,F4,uElType,uType,nodeType,AType,tType,DiffType,MinvType,StiffType})
  @femheat_deterministicpreamble
  @femheat_stochasticpreamble
  function rhs!(i,u,resid,dW,uOld)
    u = reshape(u,N,numvars)
    resid = reshape(resid,N,numvars)
    resid[freenode,:] = u[freenode,:] - uOld[freenode,:] + D.*(dt*Minv[freenode,freenode]*A[freenode,freenode]*u[freenode,:]) - (Minv*dt*quadfbasis((x,u)->f((i)*dt,x,u),(x)->gD((i)*dt,x),(x)->gN((i)*dt,x),A,u,node,elem,area,bdnode,mid,N,NT,dirichlet,neumann,islinear,numvars))[freenode,:] -(sqrtdt.*dW.*(Minv*quadfbasis((x,u)->σ((i)*dt,x,u),(x)->gD((i)*dt,x),(x)->gN((i)*dt,x),
                A,u,node,elem,area,bdnode,mid,N,NT,dirichlet,neumann,islinear,numvars)))[freenode,:]
    u = vec(u)
    resid = vec(resid)
  end
  @femheat_nonlinearsolvepreamble
  @inbounds for i=1:numiters
    t += dt
    @femheat_nonlinearsolvestochasticloop
    @femheat_footer
  end
  @femheat_postamble
end

"""
FEMmesh

Holds the information describing a finite element mesh. For information on how (node,elem)
can be interpreted as a mesh describing a geometry, see the mesh specification documentation.

### Fields

* `node`: The nodes in the (node,elem) structure.
* `elem`: The elements in the (node,elem) structure.
* `bdnode`: Vector of indices for the boundary nodes.
* `freenode`: Vector of indices for the free (non-dirichlet bound) nodes.
* `bdedge`: Indices of the edges in totaledge which are on the boundary.
* `is_bdnode`: Boolean which is true for nodes on the boundary.
* `is_bdelem`: Boolean which is true for elements on the boundary.
* `bdflag`: Flag which describes the type of boundary condition. 1=> dirichlet, 2=>neumann, 3=>robin.
* `totaledge`: Vector of the edges.
* `area`: Vector which is the area for each element.
* `dirichlet`: Indices for the nodes on the boundary which have a dirichlet boundary condition.
* `neumann`: Indices for the nodes on the boundary which have a neumann boundary condition.
* `robin`: Indices for the nodes on the boundary which have a robin boundary condition.
* `N::Int`: The number of nodes.
* `NT`::Int: The number of triangles (elements).
* `dx`: The spatial discretization size. If non-uniform, this is the average.
* `dt`: The time discretization size. If adaptive, this is the initial.
* `T`::Number: The end time.
* `numiters`::Int: The number of iterations to go from 0 to T using dt.
* `μ`: The CFL μ stability parameter.
* `ν`: The CFL ν stability parameter.
* `evolutionEq`: True for a mesh which has non-trivial time components.

"""
type FEMmesh{T1,T2,tType,TType} <: Mesh
  node::T1
  elem::Array{Int,2}
  bdnode::Vector{Int}
  freenode::Vector{Int}
  bdedge::Array{Int,2}
  is_bdnode::BitArray{1}
  is_bdelem::BitArray{1}
  bdflag::Array{Int8,2}
  totaledge::Array{Int,2}
  area::T2
  dirichlet::Array{Int,2}
  neumann::Array{Int,2}
  robin::Array{Int,2}
  N::Int
  NT::Int
  dt::tType
  T::TType
  numiters::Int
  evolutionEq::Bool
end

function FEMmesh(node,elem,dt,T,bdtype)
  N = size(node,1); NT = size(elem,1);
  totaledge = [elem[:,[2,3]]; elem[:,[3,1]]; elem[:,[1,2]]]

  #Compute the area of each element
  ve = Array{eltype(node)}(size(node[elem[:,3],:])...,3)
  ## Compute vedge, edge as a vector, and area of each element
  ve[:,:,1] = node[elem[:,3],:]-node[elem[:,2],:]
  ve[:,:,2] = node[elem[:,1],:]-node[elem[:,3],:]
  ve[:,:,3] = node[elem[:,2],:]-node[elem[:,1],:]
  area = 0.5*abs.(-ve[:,1,3].*ve[:,2,2]+ve[:,2,3].*ve[:,1,2])

  #Boundary Conditions
  bdnode,bdedge,is_bdnode,is_bdelem = findboundary(elem)
  bdflag = setboundary(node::AbstractArray,elem::AbstractArray,bdtype)
  dirichlet = totaledge[vec(bdflag .== 1),:]
  neumann = totaledge[vec(bdflag .== 2),:]
  robin = totaledge[vec(bdflag .== 3),:]
  is_bdnode = falses(N)
  is_bdnode[dirichlet] = true
  bdnode = find(is_bdnode)
  freenode = find(!is_bdnode)
  if dt != 0
    numiters = round(Int64,T/dt)
  else
    numiters = 0
  end
  FEMmesh(node,elem,bdnode,freenode,bdedge,is_bdnode,is_bdelem,bdflag,totaledge,area,dirichlet,neumann,robin,N,NT,dt,T,numiters,T!=0)
end
FEMmesh(node,elem,bdtype)=FEMmesh(node,elem,0,0,bdtype)

"""
`SimpleMesh`

Holds the information describing a finite element mesh. For information on how (node,elem)
can be interpreted as a mesh describing a geometry, see [Programming of Finite
Element Methods by Long Chen](http://www.math.uci.edu/~chenlong/226/Ch3FEMCode.pdf).

### Fields

* `node`: The nodes in the (node,elem) structure.
* `elem`: The elements in the (node,elem) structure.
"""
type SimpleMesh{T} <: Mesh
  node::T
  elem::Array{Int,2}
end


"""
`CFLμ(dt,dx)``

Computes the CFL-condition ``μ= dt/(dx*dx)``
"""
CFLμ(dt,dx)=dt/(dx*dx)

"""
`CFLν(dt,dx)``

Computes the CFL-condition ``ν= dt/dx``
"""
CFLν(dt,dx)=dt/dx

"""
`fem_squaremesh(square,h)`

Returns the grid in the iFEM form of the two arrays (node,elem)
"""
function fem_squaremesh(square,h)
  x0 = square[1]; x1= square[2];
  y0 = square[3]; y1= square[4];
  x,y = meshgrid(x0:h:x1,y0:h:y1)
  node = [x[:] y[:]];

  ni = size(x,1); # number of rows
  N = size(node,1);
  t2nidxMap = 1:N-ni;
  topNode = ni:ni:N-ni;
  t2nidxMap = deleteat!(collect(t2nidxMap),collect(topNode));
  k = t2nidxMap;
  elem = [k+ni k+ni+1 k ; k+1 k k+ni+1];
  return(node,elem)
end

"""
`notime_squaremesh(square,dx,bdtype)`

Computes the (node,elem) square mesh for the square
with the chosen `dx` and boundary settings.

###Example

```julia
square=[0 1 0 1] #Unit Square
dx=.25
notime_squaremesh(square,dx,"dirichlet")
```
"""
function notime_squaremesh(square,dx,bdtype)
  node,elem = fem_squaremesh(square,dx)
  return(FEMmesh(node,elem,bdtype))
end

"""
`parabolic_squaremesh(square,dx,dt,T,bdtype)`

Computes the `(node,elem) x [0,T]` parabolic square mesh
for the square with the chosen `dx` and boundary settings
and with the constant time intervals `dt`.

###Example

```julia
square=[0 1 0 1] #Unit Square
dx=.25; dt=.25;T=2
parabolic_squaremesh(square,dx,dt,T,:dirichlet)
```

"""
function parabolic_squaremesh(square,dx,dt,T,bdtype)
  node,elem = fem_squaremesh(square,dx)
  return(FEMmesh(node,elem,dt,T,bdtype))
end

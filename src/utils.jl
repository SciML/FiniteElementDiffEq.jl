"""
`getNoise(N,node,elem;noisetype=:White)`

Returns a random vector corresponding to the noise type which was chosen.
"""
function getNoise(u,node,elem;noisetype=:White)
  if noisetype==:White
    return(ChunkedArray(randn,u))
  end
end

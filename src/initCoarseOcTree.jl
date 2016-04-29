export initCoarseOcTree

function initCoarseOcTree(Lims::Array{Float64,2},nfine::Vector{Int64},coarseLev::Int64)
  #initCoarseOcTree(Lims,nfine,ncoarse)
  #Creates a uniform OcTree with 
  #2^(nfine[1]-coarseLev) X 2^(nfine[2]-coarseLev) X 2^(nfine[3]-coarseLev) cells on
  #an underlying mesh of 2^nfine[1] X 2^nfine[2] X 2^nfine[3] fine cells
  
  @assert all(nfine .>= 1) "initCoarseOcTree: fine mesh must have at least 2 cells in all dimensions"
  @assert size(Lims) == (2,3) "initCoarseOcTree: Lims must have size (2,3)"
  nCoarse = nfine - coarseLev
  @assert all(nCoarse .>= 0) "Coarsening level too high for underlying mesh in at least 1 dimension."
  
  Lx = Lims[2,1] - Lims[1,1]
  Ly = Lims[2,2] - Lims[1,2]
  Lz = Lims[2,3] - Lims[1,3]
  h  = [Lx;Ly;Lz]./(2.^nfine)
  
  i,j,k = ndgrid(1:2^coarseLev:2^nfine[1],1:2^coarseLev:2^nfine[2],1:2^coarseLev:2^nfine[3])
  ncls  = prod(size(i))
  bsz   = (2^coarseLev)*ones(Int64,ncls)
  S     = sparse3(vec(i),vec(j),vec(k),bsz,2.^nfine)
  
  return S,h
end
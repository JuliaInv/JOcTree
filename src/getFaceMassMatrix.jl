export getFaceMassMatrix, getdFaceMassMatrix

function getFaceMassMatrix(M::OcTreeMeshFV,sigma::Vector)
  # For octree meshes
  if isempty(M.Pf)
    M.Pf = getFaceMassMatrixIntegrationMatrix(M.S, M.h)
  end
  P = M.Pf
  n = length(sigma)
  if n == M.nc
    M = P' * kron(speye(24),spdiagm(sigma)) * P
  elseif n == 3 * M.nc
    M = P' * kron(speye(8),spdiagm(sigma)) * P
  elseif n == 6 * M.nc
    R12 = kron(sparse([2,1,3],[1,2,3],1),speye(Int,M.nc))
    R23 = kron(sparse([1,3,2],[1,2,3],1),speye(Int,M.nc))
    D   = spdiagm(sigma[       1:3*M.nc])
    N   = spdiagm(sigma[3*M.nc+1:6*M.nc])
    S   = D + R12 * N * R23 + R23 * N * R12
    M   = P' * kron(speye(8),S) * P
   else
     error("Invalid size")
   end
	 return M
end

# For backwards compatibility
getdFaceMassMatrix(M::OcTreeMeshFV,v::Vector) = getdFaceMassMatrix(M,ones(M.nc),v)

function getdFaceMassMatrix(M::OcTreeMeshFV,sigma::Vector,v::Vector)
  # Derivative
  if isempty(M.Pf)
    M.Pf = getFaceMassMatrixIntegrationMatrix(M.S, M.h)
  end
  P = M.Pf
  w = P * v
  n = length(sigma)
  if n == M.nc
    dM = P' * spdiagm(w) * kron(ones(24),speye(M.nc))
  elseif n == 3 * M.nc
    dM = P' * spdiagm(w) * kron(ones(8),speye(3*M.nc))
  elseif n == 6 * M.nc
    R12 = kron(speye(Int,8),kron(sparse([2,1,3],[1,2,3],1),speye(Int,M.nc)))
    R23 = kron(speye(Int,8),kron(sparse([1,3,2],[1,2,3],1),speye(Int,M.nc)))
    D   = spdiagm(w)
    N   = R12 * spdiagm(R23 * w) + R23 * spdiagm(R12 * w)
    dM  = hcat(
      P' * D * kron(ones(8),speye(3*M.nc)),
      P' * N * kron(ones(8),speye(3*M.nc)))
  else
    error("Invalid size")
  end
	return dM
end

function getFaceMassMatrixIntegrationMatrix(S::SparseArray3D,h)
    
    n   = S.sz;
    nex = (n[1]+1,n[2],n[3])
    ney = (n[1],n[2]+1,n[3])
    nez = (n[1],n[2],n[3]+1)
    
    i,j,k,bsz = find3(S)
    ex,ey,ez = getFaceNumbering(S)
    
    nx = nnz(ex)
    ny = nnz(ey)
    nz = nnz(ez)
    nc = nnz(S)
    c  = collect(1:nc)
    uc = sqrt(bsz.^3 * prod(h) / 8)
    
    Px1 = getNodesFromIndices(ex.SV,nex,i,j,k)
    Py1 = getNodesFromIndices(ey.SV,ney,i,j,k)
    Pz1 = getNodesFromIndices(ez.SV,nez,i,j,k)
    
    Px2 = getNodesFromIndices(ex.SV,nex,i+bsz,j,k)
    Py2 = getNodesFromIndices(ey.SV,ney,i,j,k)
    Pz2 = getNodesFromIndices(ez.SV,nez,i,j,k)
    
    Px3 = getNodesFromIndices(ex.SV,nex,i,j,k)
    Py3 = getNodesFromIndices(ey.SV,ney,i,j+bsz,k)
    Pz3 = getNodesFromIndices(ez.SV,nez,i,j,k)
    
    Px4 = getNodesFromIndices(ex.SV,nex,i+bsz,j,k)
    Py4 = getNodesFromIndices(ey.SV,ney,i,j+bsz,k)
    Pz4 = getNodesFromIndices(ez.SV,nez,i,j,k)
    
    Px5 = getNodesFromIndices(ex.SV,nex,i,j,k)
    Py5 = getNodesFromIndices(ey.SV,ney,i,j,k)
    Pz5 = getNodesFromIndices(ez.SV,nez,i,j,k+bsz)
    
    Px6 = getNodesFromIndices(ex.SV,nex,i+bsz,j,k)
    Py6 = getNodesFromIndices(ey.SV,ney,i,j,k)
    Pz6 = getNodesFromIndices(ez.SV,nez,i,j,k+bsz)
    
    Px7 = getNodesFromIndices(ex.SV,nex,i,j,k)
    Py7 = getNodesFromIndices(ey.SV,ney,i,j+bsz,k)
    Pz7 = getNodesFromIndices(ez.SV,nez,i,j,k+bsz)
    
    Px8 = getNodesFromIndices(ex.SV,nex,i+bsz,j,k)
    Py8 = getNodesFromIndices(ey.SV,ney,i,j+bsz,k)
    Pz8 = getNodesFromIndices(ez.SV,nez,i,j,k+bsz)
    
    sp1(j) = sparse(c,j,uc,nc,nx)
    sp2(j) = sparse(c,j,uc,nc,ny)
    sp3(j) = sparse(c,j,uc,nc,nz)
    
    P1 = blkdiag(sp1(Px1),sp2(Py1),sp3(Pz1))
    P2 = blkdiag(sp1(Px2),sp2(Py2),sp3(Pz2))
    P3 = blkdiag(sp1(Px3),sp2(Py3),sp3(Pz3))
    P4 = blkdiag(sp1(Px4),sp2(Py4),sp3(Pz4))
    P5 = blkdiag(sp1(Px5),sp2(Py5),sp3(Pz5))
    P6 = blkdiag(sp1(Px6),sp2(Py6),sp3(Pz6))
    P7 = blkdiag(sp1(Px7),sp2(Py7),sp3(Pz7))
    P8 = blkdiag(sp1(Px8),sp2(Py8),sp3(Pz8))
    
    P = vcat(P1, P2, P3, P4, P5, P6, P7, P8)
    
    return P
    
end

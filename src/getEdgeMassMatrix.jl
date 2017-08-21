
export getEdgeMassMatrix, getdEdgeMassMatrix,
       dEdgeMassMatrixTimesVector, dEdgeMassMatrixTrTimesVector


"""
function getEdgeMassMatrix(M::OcTreeMeshFV,sigma::Vector)

Returns finite volume edge mass matrix. Can handle isotropic and anisotropic models.


Input:

   M::OcTreeMeshFV Octree mesh object
   sigma::vector Conductivity model. Can be used for isotropic and anisotropic models depending on length
          of sigma. For an isotropic model length(sigma) = nc, where nc is the number of mesh cells.
          Use length(sigma)=3*nc a diagonally anisotropic model and length(sigma)=6*nc for
          for fully anisotropic model.
Output:
   Me::SparseMatrixCSC Edge mass matrix

"""
function getEdgeMassMatrix(Mesh::OcTreeMeshFV,sigma::Vector)
    # For octree meshes

    n = length(sigma)
    @assert in(n,[Mesh.nc;3*Mesh.nc;6*Mesh.nc]) "Invalid size of sigma"
    if n == Mesh.nc
        #M = Pt * kron(speye(24),spdiagm(sigma)) * P
        Ae,Aet = getEdgeAverageMatrix(Mesh)
        v      = getVolumeVector(Mesh)
        m      = Vector{eltype(sigma)}(sum(Mesh.ne))
        A_mul_B!(m,Aet,3*v.*sigma)
        M  = spdiagm(m)  #M = spdiagm(Ae'*(v.*sigma))
    else
        if isempty(Mesh.Pe)
            Mesh.Pe  = getEdgeMassMatrixIntegrationMatrix(Mesh.S, Mesh.h)
            Mesh.Pet = Mesh.Pe'
        end
        P  = Mesh.Pe
        Pt = Mesh.Pet
    end
    if n == 3 * Mesh.nc
        M = Pt * kron(speye(8),spdiagm(sigma)) * P
    elseif n == 6 * Mesh.nc
        R12 = kron(sparse([2,1,3],[1,2,3],1),speye(Int,Mesh.nc))
        R23 = kron(sparse([1,3,2],[1,2,3],1),speye(Int,Mesh.nc))
        D   = spdiagm(sigma[       1:3*Mesh.nc])
        N   = Diagonal(sigma[3*Mesh.nc+1:6*Mesh.nc])
        S   = D + R12 * N * R23 + R23 * N * R12
        M   = Pt * kron(speye(8),S) * P
    end
    return M
end

function getdEdgeMassMatrix{T<:Number}(Mesh::OcTreeMeshFV, sigma::Vector, v::Vector{T})
    # Derivative

    n = length(sigma)
    @assert in(n,[Mesh.nc;3*Mesh.nc;6*Mesh.nc]) "Invalid size of sigma"
    if n == Mesh.nc
        Ae,Aet = getEdgeAverageMatrix(Mesh)
        vol    = getVolumeVector(Mesh)
        dM = SparseMatrixCSC(Aet.m, Aet.n, copy(Aet.colptr),
                             copy(Aet.rowval), T.(3.*Aet.nzval))
        DiagTimesMTimesDiag!(v,dM,vol)
    else
        if isempty(Mesh.Pe)
            Mesh.Pe  = getEdgeMassMatrixIntegrationMatrix(Mesh.S, Mesh.h)
            Mesh.Pet = Mesh.Pe'
        end
        P  = Mesh.Pe
        Pt = Mesh.Pet
        w = P * v
    end
    if n == 3 * Mesh.nc
        K  = kron(ones(T, 8),speye(3*Mesh.nc))
        dM = Pt * DiagTimesM(w,K)
    elseif n == 6 * Mesh.nc
        R12 = kron(speye(Int,8),kron(sparse([2,1,3],[1,2,3],1),speye(Int,Mesh.nc)))
        R23 = kron(speye(Int,8),kron(sparse([1,3,2],[1,2,3],1),speye(Int,Mesh.nc)))
        D   = spdiagm(w)
        N   = R12 * spdiagm(R23 * w) + R23 * spdiagm(R12 * w)
        dM  = hcat(
        Pt * D * kron(ones(8),speye(3*Mesh.nc)),
        Pt * N * kron(ones(8),speye(3*Mesh.nc)))
    end
	return dM
end

# for backwards compatibility
getdEdgeMassMatrix(M::OcTreeMeshFV,v::Vector) = getdEdgeMassMatrix(M, ones(M.nc), v)

function getEdgeMassMatrixIntegrationMatrix(S, h)
# P = getEdgeMassMatrixIntegrationMatrix(S, h)
# Compute a matrix P such that the edge function mass matrix
# for a general symmetric anisotropic coefficient can be computed by
#   M = P' * C * P
# Here,
#       | diag(cxx) diag(cxy) diag(cxz) |
#   C = | diag(cxy) diag(cyy) diag(cyz) |
#       | diag(cxz) diag(cyz) diag(czz) |
# expands the six components cxx, cyy, czz, cxy, cxz, cyz of the symmetric
# coefficient tensor for all numCell cells to a [3*numCells]^2 matrix.
#
# To illustrate the integration, consider the 2D case (QuadTree) and, in
# particular, the cell C1 of size 2 * h which has two right neighbors of
# of size h:
#
#            ex2
#      +-------------+------+------+
#      |q3         q4|      |      |
#      |             | ey3  |      |
#      |             |      |      |
#  ey1 |      C1     +------+------+
#      |             |      |      |
#      |             | ey2  |      |
#      |q1         q2|      |      |
#      +-------------+------+------+
#            ex1
#
# The integration of the mass matrix for the cell is performed in three
# steps:
# 1. Interpolate the edge function, which contains the x- and y-field
#    components ex1, ex2, ey1, ey2, ey3 staggered on the x- and y- edges,
#    to the four quadrature points qi (i = 1...4)
#      Ex(qi) = Px(ex1,ex2)
#      Ey(qi) = Py(ey1,ey2,ey3)
# 2. Compute the inner product of the collocated field, weighted by the
#    tensor C1
#      I1 = (Ex(qi), Ey(qi)) * C1 * (Ex(qi), Ey(qi))^T   (i = 1...4)
# 3. Integrate using numerical quadrature
#      M1 = (I1 + I2 + I3 + I4) / 4 * (2 * h)^2
#
# We use the following nearest neighbour interpolation
#   Ex(q1) = ex1, Ey(q1) = ey1
#   Ex(q2) = ex1, Ey(q2) = ey2
#   Ex(q3) = ex2, Ey(q3) = ey1
#   Ex(q4) = ex2, Ey(q4) = ey3
#
# The mass matrix has up to 9 non-zero entries per row/column.
#
# Note that this integration method neglects hanging edges. The resulting
# mass matrix M for an isotropic medium differs from the mass matrix Me
# which is obtained by the standard integration method by
#   Ae = getEdgeToCellCenterMatrix(S); # edge to cell center average
#   v  = prod(h) * nonzeros(S).^3;     # cell volume
#   Me = sdiag( (v .* c)' * Ae );
# where c is the isotropic coefficient for every cell. Let
#   C  = blkdiag(sdiag(c), sdiag(c), sdiag(c));
#   M  = P' * C * P;
# Then,
#   norm(M - Me) > 0
# unless the OcTree is degenerate (a regular mesh).

# C. Schwarzbach, April 2014

# edge numbering (implies location of edge)

EX,EY,EZ  = getEdgeNumbering(S)

i,j,k,bsz = find3(S)
mx,my,mz  = size(S)
nex = nnz(EX)
ney = nnz(EY)
nez = nnz(EZ)

# locate edge numbers for quadrature points

# x-edges
i0 = i
i1 = floor.(Integer,i + bsz / 2)
j0 = j
j1 = j + bsz
k0 = k
k1 = k + bsz

ex000 = getNodesFromIndices(EX.SV,(mx,my+1,mz+1),i0,j0,k0)
ex100 = getNodesFromIndices(EX.SV,(mx,my+1,mz+1),i1,j0,k0)
ex010 = getNodesFromIndices(EX.SV,(mx,my+1,mz+1),i0,j1,k0)
ex110 = getNodesFromIndices(EX.SV,(mx,my+1,mz+1),i1,j1,k0)
ex001 = getNodesFromIndices(EX.SV,(mx,my+1,mz+1),i0,j0,k1)
ex101 = getNodesFromIndices(EX.SV,(mx,my+1,mz+1),i1,j0,k1)
ex011 = getNodesFromIndices(EX.SV,(mx,my+1,mz+1),i0,j1,k1)
ex111 = getNodesFromIndices(EX.SV,(mx,my+1,mz+1),i1,j1,k1)

ex100 = merge(ex100, ex000)
ex110 = merge(ex110, ex010)
ex101 = merge(ex101, ex001)
ex111 = merge(ex111, ex011)

# y-edges
i0 = i
i1 = i + bsz
j0 = j
j1 = floor.(Integer,j + bsz / 2)
k0 = k
k1 = k + bsz

ey000 = getNodesFromIndices(EY.SV,(mx+1,my,mz+1),i0,j0,k0)
ey100 = getNodesFromIndices(EY.SV,(mx+1,my,mz+1),i1,j0,k0)
ey010 = getNodesFromIndices(EY.SV,(mx+1,my,mz+1),i0,j1,k0)
ey110 = getNodesFromIndices(EY.SV,(mx+1,my,mz+1),i1,j1,k0)
ey001 = getNodesFromIndices(EY.SV,(mx+1,my,mz+1),i0,j0,k1)
ey101 = getNodesFromIndices(EY.SV,(mx+1,my,mz+1),i1,j0,k1)
ey011 = getNodesFromIndices(EY.SV,(mx+1,my,mz+1),i0,j1,k1)
ey111 = getNodesFromIndices(EY.SV,(mx+1,my,mz+1),i1,j1,k1)

ey010 = merge(ey010, ey000)
ey110 = merge(ey110, ey100)
ey011 = merge(ey011, ey001)
ey111 = merge(ey111, ey101)

# z-edges
i0 = i
i1 = i + bsz
j0 = j
j1 = j + bsz
k0 = k
k1 = floor.(Integer,k + bsz / 2)

ez000 = getNodesFromIndices(EZ.SV,(mx+1,my+1,mz),i0,j0,k0)
ez100 = getNodesFromIndices(EZ.SV,(mx+1,my+1,mz),i1,j0,k0)
ez010 = getNodesFromIndices(EZ.SV,(mx+1,my+1,mz),i0,j1,k0)
ez110 = getNodesFromIndices(EZ.SV,(mx+1,my+1,mz),i1,j1,k0)
ez001 = getNodesFromIndices(EZ.SV,(mx+1,my+1,mz),i0,j0,k1)
ez101 = getNodesFromIndices(EZ.SV,(mx+1,my+1,mz),i1,j0,k1)
ez011 = getNodesFromIndices(EZ.SV,(mx+1,my+1,mz),i0,j1,k1)
ez111 = getNodesFromIndices(EZ.SV,(mx+1,my+1,mz),i1,j1,k1)

ez001 = merge(ez001, ez000)
ez101 = merge(ez101, ez100)
ez011 = merge(ez011, ez010)
ez111 = merge(ez111, ez110)

# cell indices
nc = length(bsz)
c  = collect(1:nc)

# set values for each cell to square root of cell volume and quadrature
# weight; product of P' * C * P gives the correct scaling by cell volume
# and quadrature weight
uc = sqrt.(bsz.^3 * prod(h) / 8)

# integrate using edge to cell corner interpolation
Px = sparse(c, ex000, uc, nc, nex) # Px1
Py = sparse(c, ey000, uc, nc, ney) # Py1
Pz = sparse(c, ez000, uc, nc, nez) # Pz1
P1 = blkdiag(Px, Py, Pz)

Px = sparse(c, ex100, uc, nc, nex) # Px2
Py = sparse(c, ey100, uc, nc, ney) # Py2
Pz = sparse(c, ez100, uc, nc, nez) # Pz2
P2 = blkdiag(Px, Py, Pz)

Px = sparse(c, ex010, uc, nc, nex) # Px3
Py = sparse(c, ey010, uc, nc, ney) # Py3
Pz = sparse(c, ez010, uc, nc, nez) # Pz3
P3 = blkdiag(Px, Py, Pz)

Px = sparse(c, ex110, uc, nc, nex) # Px4
Py = sparse(c, ey110, uc, nc, ney) # Py4
Pz = sparse(c, ez110, uc, nc, nez) # Pz4
P4 = blkdiag(Px, Py, Pz)

Px = sparse(c, ex001, uc, nc, nex) # Px5
Py = sparse(c, ey001, uc, nc, ney) # Py5
Pz = sparse(c, ez001, uc, nc, nez) # Pz5
P5 = blkdiag(Px, Py, Pz)

Px = sparse(c, ex101, uc, nc, nex) # Px6
Py = sparse(c, ey101, uc, nc, ney) # Py6
Pz = sparse(c, ez101, uc, nc, nez) # Pz6
P6 = blkdiag(Px, Py, Pz)

Px = sparse(c, ex011, uc, nc, nex) # Px7
Py = sparse(c, ey011, uc, nc, ney) # Py7
Pz = sparse(c, ez011, uc, nc, nez) # Pz7
P7 = blkdiag(Px, Py, Pz)

Px = sparse(c, ex111, uc, nc, nex) # Px8
Py = sparse(c, ey111, uc, nc, ney) # Py8
Pz = sparse(c, ez111, uc, nc, nez) # Pz8
P8 = blkdiag(Px, Py, Pz)

# Collect everything
P = vcat(P1, P2, P3, P4, P5, P6, P7, P8)

return P

end

function merge(a, b)
# copy entries of b into a for zero entries of a
c      = copy(a)
idz    = c .== 0
c[idz] = b[idz]
return c
end

function dEdgeMassMatrixTimesVector(Mesh::OcTreeMeshFV, sigma::Vector, v::Vector, x::Vector)
    # Derivative (getdEdgeMassMatrix) times a vector(x)

    n = length(sigma)
    if length(x) != n
        error("length(x) != length(sigma)")
    end
    @assert in(n,[Mesh.nc;3*Mesh.nc;6*Mesh.nc]) "Invalid size of sigma"

    if n > Mesh.nc
        if isempty(Mesh.Pe)
            Mesh.Pe = getEdgeMassMatrixIntegrationMatrix(Mesh.S, Mesh.h)
        end
    end

    if n == 6 * Mesh.nc
        # Not the best solution!
        dM = getdEdgeMassMatrix(Mesh, sigma, v)
        return dM * x
    end

    if n == Mesh.nc
        vol = getVolumeVector(Mesh)
        Ae,Aet = getEdgeAverageMatrix(Mesh)
        volx = vol.*x
        Aetvolx = 3*Aet*volx
        dMx = v.*Aetvolx
    elseif n == 3 * Mesh.nc
        pv  = Mesh.Pe * v
        dMx = Mesh.Pe' * (pv .* repmat(x, 8))
    end

    return dMx
end


function dEdgeMassMatrixTrTimesVector(Mesh::OcTreeMeshFV, sigma::Vector, v::Vector, x::Vector)
    # Derivative (getdEdgeMassMatrix) transpose times a vector(x)


    n = length(sigma)
    @assert in(n,[Mesh.nc;3*Mesh.nc;6*Mesh.nc]) "Invalid size of sigma"
    if n > Mesh.nc
        if isempty(Mesh.Pe)
            Mesh.Pe = getEdgeMassMatrixIntegrationMatrix(Mesh.S,Mesh.h)
        end
    end

    if n == 6 * Mesh.nc
        # Not the best solution!
        dM = getdEdgeMassMatrix(Mesh, sigma, v)
        return dM' * x
    end

    # dM' = kron(ones(24,1),speye(nnz(M.S)))' * sdiag(M.Pe*v) * M.Pe

    if n == Mesh.nc
        vol  = getVolumeVector(Mesh)
        Ae,  = getEdgeAverageMatrix(Mesh)
        vx   = v.*x
        Aevx = 3*Ae*vx
        dMTx = vol.*Aevx
    elseif n == 3 * Mesh.nc
        pv   = Mesh.Pe * v
        dd   = pv .* (Mesh.Pe * conj(x))
        pv   = []
        ns   = 3 * Mesh.nc
        i2   = 8
        dMTx = dd[1:ns]
        j    = ns
        for i = 2:i2
            @inbounds dMTx += dd[j+1 : j+ns]
            j += ns
        end  # i
        dMTx = conj(dMTx)
    end





    return dMTx
end

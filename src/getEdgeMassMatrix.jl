export getEdgeMassMatrix, getdEdgeMassMatrix,
       dEdgeMassMatrixTimesVector, dEdgeMassMatrixTrTimesVector

"""
    M = getEdgeMassMatrix(mesh::OcTreeMeshFV, sigma::Vector)

Compute finite volume edge mass matrix with coefficient `sigma` on OcTree
mesh `mesh`.

`getEdgeMassMatrix` distinguishes three cases:
 - `length(sigma) == mesh.nc`: scalar coefficient (isotropy)
 - `length(sigma) == mesh.nc * 3`: diagonal tensor coefficient (diagonal anisotropy)
 - `length(sigma) == mesh.nc * 6`: symmetric tensor coefficient (general anisotropy)

The tensor coefficients must be stored in the following order:
 - `sigma = vcat(sigma_xx, sigma_yy, sigma_zz)`
 - `sigma = vcat(sigma_xx, sigma_yy, sigma_zz, sigma_xy, sigma_xz, sigma_yz)`
where the vectors `sigma_xx`, ... contain the respective entries of the
tensor for all OcTree cells.

The first time that one of `getEdgeMassMatrix`, `getdEdgeMassMatrix`,
`dEdgeMassMatrixTimesVector` or `dEdgeMassMatrixTrTimesVector` is called for
the mesh and for a particular type of coefficient, integration weights and
data structure are precomputed and stored internally. Subsequent calls for
the same type of coefficient execute much faster by exploiting the stored
data structure.

The mass matrix integration neglects hanging edges. For best results,
elimiate hanging edges as follows:

    M = getEdgeMassMatrix(mesh, sigma)
    N, = getEdgeConstraints(mesh)
    M = N' * M * N
"""
function getEdgeMassMatrix(mesh::OcTreeMeshFV,sigma::Vector)
    P = setupEdgeMassMatrix(mesh, sigma)
    nzval = P.A * sigma
    M = SparseMatrixCSC(P.n, P.n, P.colptr, P.rowval, nzval)
    return M
end

"""
    dM = getdEdgeMassMatrix(mesh::OcTreeMeshFV, sigma::Vector, v::Vector)

Compute directional derivative of edge mass matrix w.r.t. coefficient `sigma`
```math
\\nabla_\\sigma (M(\\sigma) v)
```

`length(v)` must equal the number of edges.

`getdEdgeMassMatrix` returns of sparse matrix of size
`(length(v), length(sigma))`.

See documentation of `getEdgeMassMatrix` for more details.
"""
function getdEdgeMassMatrix{T<:Number}(mesh::OcTreeMeshFV, sigma::Vector, v::Vector{T})
    P = setupEdgeMassMatrix(mesh, sigma)
    nz = length(P.colval)
    if nz == 0 # isotropic or diagonally anisotropic
        dM = SparseMatrixCSC(P.A.m, P.A.n, copy(P.A.colptr), copy(P.A.rowval), T.(P.A.nzval))
        DiagTimesM(v, dM)
    else # generally anisotropic
        N = eltype(P.A.colptr)
        V = SparseMatrixCSC(P.n, nz, collect(N,1:nz+1), P.rowval, v[P.colval])
        dM = V * P.A
    end
    return dM
end

# for backwards compatibility
getdEdgeMassMatrix(mesh::OcTreeMeshFV, v::Vector) = getdEdgeMassMatrix(mesh, ones(mesh.nc), v)

"""
    dMx = dEdgeMassMatrixTimesVector(mesh::OcTreeMeshFV, sigma::Vector, v::Vector, x::Vector)

Compute matrix-vector product between directional derivative of edge mass matrix
 w.r.t. coefficient `sigma`
```math
\\nabla_\\sigma (M(\\sigma) v)
```
and vector `x`.

`dEdgeMassMatrixTimesVector(mesh, sigma, v, x)` is a fast implementation of
`getdEdgeMassMatrix(mesh, sigma, v) * x`.

`length(v)` must equal the number of edges;
`x` must satisfy `length(x) == length(sigma)`.

`dEdgeMassMatrixTimesVector` returns a vector of size `length(v)`.

See documentation of `getEdgeMassMatrix` for more details.
"""
function dEdgeMassMatrixTimesVector{T<:Number}(mesh::OcTreeMeshFV, sigma::Vector, v::Vector{T}, x::Vector)
    P = setupEdgeMassMatrix(mesh, sigma)
    nz = length(P.colval)
    if nz == 0 # isotropic or diagonally anisotropic
        dMx = T.(P.A * x)
        dMx .*= v
    else # generally anisotropic
        V = SparseMatrixCSC(P.n, nz, collect(eltype(P.rowval),1:nz+1), P.rowval, v[P.colval])
        dMx = V * (P.A * x)
    end
    return dMx
end

"""
    dMTx = dEdgeMassMatrixTrTimesVector(mesh::OcTreeMeshFV, sigma::Vector, v::Vector, x::Vector)

Compute matrix-vector product between (complex) transpose of  directional
derivative of edge mass matrix
 w.r.t. coefficient `sigma`
```math
\\nabla_\\sigma (M(\\sigma) v)
```
and vector `x`.

`dEdgeMassMatrixTimesTrVector(mesh, sigma, v, x)` is a fast implementation of
`getdEdgeMassMatrix(mesh, sigma, v)' * x`.

`length(v)` must equal the number of edges;
`x` must satisfy `length(x) == length(v)`.

`dEdgeMassMatrixTimesTrVector` returns a vector of size `length(sigma)`.

See documentation of `getEdgeMassMatrix` for more details.
"""
function dEdgeMassMatrixTrTimesVector{T<:Number}(mesh::OcTreeMeshFV, sigma::Vector, v::Vector{T}, x::Vector{T})
    P = setupEdgeMassMatrix(mesh, sigma)
    nz = length(P.colval)
    if nz == 0 # isotropic or diagonally anisotropic
        u = conj.(v) .* x
        dMTx = P.A' * u
    else # generally anisotropic
        V = SparseMatrixCSC(P.n, nz, collect(eltype(P.rowval),1:nz+1), P.rowval, v[P.colval])
        u = V' * x
        dMTx = P.A' * u
    end
    return dMTx
end

### Private (not exported) methods below

"""
    P = setupEdgeMassMatrix(mesh, sigma)

Precompute and return stored data structure for edge mass matrix integration
for isotropic, diagonally anisotropic or generally anisotropic coefficient.
"""
function setupEdgeMassMatrix(mesh::OcTreeMeshFV, sigma)

    na = length(sigma)
    nc = mesh.nc
    N  = typeof(nc)

    @assert in(na, [nc, 3*nc, 6*nc]) "Invalid size of sigma"

    if !haskey(mesh.Pe, na)

        n = sum(mesh.ne)
        ex, ey, ez, w = getEdgeMassMatrixQuadrature(mesh)

        if na == nc

            i = vcat(ex, ey, ez)
            j = repmat(one(N):nc, 24)
            a = repmat(w, 24)
            A = sparse(i, j, a, n, na)
            colptr = collect(N,1:n+1)
            rowval = collect(N,1:n)
            colval = Array{N}(0) # unused

        elseif na == 3 * nc

            i = vcat(ex, ey, ez)
            j = vcat(repmat(UnitRange{N}(     1:nc  ), 8),
                     repmat(UnitRange{N}(  nc+1:2*nc), 8),
                     repmat(UnitRange{N}(2*nc+1:3*nc), 8))
            a = repmat(w, 24)
            A = sparse(i, j, a, n, na)
            colptr = collect(N,1:n+1)
            rowval = collect(N,1:n)
            colval = Array{N}(0) # unused

        elseif na == 6 * nc

            # find nonzero pattern of mass matrix
            i = vcat(ex, ey, ez, ex, ex, ey, ey, ez, ez)
            j = vcat(ex, ey, ez, ey, ez, ez, ex, ex, ey)
            P = sparse(i, j, j, n, n, (x,y) -> x)
            colptr = P.colptr
            rowval = P.rowval
            colval = copy(P.nzval) # P.nzval will be overwritten in next step

            # find destination of each quadrature point in sparse
            # mass matrix vector of nonzeros; overwrite array i
            nz = length(P.nzval)
            P.nzval .= 1:nz
            for m = 1:length(i)
              i[m] = P[i[m],j[m]]
            end

            # construct mapping of anisotropic cell property to
            # nonzero entries in mass matrix
            j = vcat(repmat(UnitRange{N}(     1:  nc), 8),
                     repmat(UnitRange{N}(  nc+1:2*nc), 8),
                     repmat(UnitRange{N}(2*nc+1:3*nc), 8),
                     repmat(UnitRange{N}(3*nc+1:4*nc), 8),
                     repmat(UnitRange{N}(4*nc+1:5*nc), 8),
                     repmat(UnitRange{N}(5*nc+1:6*nc), 8),
                     repmat(UnitRange{N}(3*nc+1:4*nc), 8),
                     repmat(UnitRange{N}(4*nc+1:5*nc), 8),
                     repmat(UnitRange{N}(5*nc+1:6*nc), 8))
            a = repmat(w, 72)
            A = sparse(i, j, a, nz, na)

        end

        mesh.Pe[na] = MassMatrix(n, A, rowval, colptr, colval)

    end

    return mesh.Pe[na]

end

"""
    ex, ey, ez, w = getEdgeMassMatrixQuadrature(mesh)

Compute quadrature points and weights for integration of edge mass matrix.

`getEdgeMassMatrixQuadrature` returns a list of edges and weights
for eight quadrature points, the corners of each cell.
"""
function getEdgeMassMatrixQuadrature(mesh::OcTreeMeshFV)

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
# Note that this integration method neglects hanging edges. The
# mass matrix
#   M1 = getEdgeMassMatrix(mesh, sigma)
# for isotropic sigma differs from
#   Ae = getEdgeToCellCenterMatrix(mesh) # edge to cell center average
#   V  = getVolume(mesh)                 # cell volume
#   M2 = sdiag(Ae' * (V * sigma))
# unless the OcTree degenerates to a regular mesh.
#
# C. Schwarzbach, April 2014

# edge numbering (implies location of edge)
EX,EY,EZ  = getEdgeNumbering(mesh)
nex = nnz(EX)
ney = nnz(EY)
nez = nnz(EZ)

# cells
i,j,k,bsz = find3(mesh.S)
mx,my,mz  = size(mesh.S)

# locate edge numbers for quadrature points
bsz2 = div.(bsz, 2)
i0 = i
i1 = i + bsz
i2 = i + bsz2
j0 = j
j1 = j + bsz
j2 = j + bsz2
k0 = k
k1 = k + bsz
k2 = k + bsz2

# x-edges
ex000 = getNodesFromIndices(EX.SV,(mx,my+1,mz+1),i0,j0,k0)
ex100 = getNodesFromIndices(EX.SV,(mx,my+1,mz+1),i2,j0,k0)
ex010 = getNodesFromIndices(EX.SV,(mx,my+1,mz+1),i0,j1,k0)
ex110 = getNodesFromIndices(EX.SV,(mx,my+1,mz+1),i2,j1,k0)
ex001 = getNodesFromIndices(EX.SV,(mx,my+1,mz+1),i0,j0,k1)
ex101 = getNodesFromIndices(EX.SV,(mx,my+1,mz+1),i2,j0,k1)
ex011 = getNodesFromIndices(EX.SV,(mx,my+1,mz+1),i0,j1,k1)
ex111 = getNodesFromIndices(EX.SV,(mx,my+1,mz+1),i2,j1,k1)

merge!(ex100, ex000)
merge!(ex110, ex010)
merge!(ex101, ex001)
merge!(ex111, ex011)

# y-edges
ey000 = getNodesFromIndices(EY.SV,(mx+1,my,mz+1),i0,j0,k0)
ey100 = getNodesFromIndices(EY.SV,(mx+1,my,mz+1),i1,j0,k0)
ey010 = getNodesFromIndices(EY.SV,(mx+1,my,mz+1),i0,j2,k0)
ey110 = getNodesFromIndices(EY.SV,(mx+1,my,mz+1),i1,j2,k0)
ey001 = getNodesFromIndices(EY.SV,(mx+1,my,mz+1),i0,j0,k1)
ey101 = getNodesFromIndices(EY.SV,(mx+1,my,mz+1),i1,j0,k1)
ey011 = getNodesFromIndices(EY.SV,(mx+1,my,mz+1),i0,j2,k1)
ey111 = getNodesFromIndices(EY.SV,(mx+1,my,mz+1),i1,j2,k1)

merge!(ey010, ey000)
merge!(ey110, ey100)
merge!(ey011, ey001)
merge!(ey111, ey101)

# z-edges
ez000 = getNodesFromIndices(EZ.SV,(mx+1,my+1,mz),i0,j0,k0)
ez100 = getNodesFromIndices(EZ.SV,(mx+1,my+1,mz),i1,j0,k0)
ez010 = getNodesFromIndices(EZ.SV,(mx+1,my+1,mz),i0,j1,k0)
ez110 = getNodesFromIndices(EZ.SV,(mx+1,my+1,mz),i1,j1,k0)
ez001 = getNodesFromIndices(EZ.SV,(mx+1,my+1,mz),i0,j0,k2)
ez101 = getNodesFromIndices(EZ.SV,(mx+1,my+1,mz),i1,j0,k2)
ez011 = getNodesFromIndices(EZ.SV,(mx+1,my+1,mz),i0,j1,k2)
ez111 = getNodesFromIndices(EZ.SV,(mx+1,my+1,mz),i1,j1,k2)

merge!(ez001, ez000)
merge!(ez101, ez100)
merge!(ez011, ez010)
merge!(ez111, ez110)

# collect everything
ex = vcat(ex000,ex100,ex010,ex110,ex001,ex101,ex011,ex111)
ey = vcat(ey000,ey100,ey010,ey110,ey001,ey101,ey011,ey111)
ez = vcat(ez000,ez100,ez010,ez110,ez001,ez101,ez011,ez111)

# shift ey and ez
ey .+= nex
ez .+= (nex + ney)

# integration weights
w = bsz.^3 * (prod(mesh.h) / 8)

return ex, ey, ez, w

end

export getNodalMassMatrix, getdNodalMassMatrix,
       dNodalMassMatrixTimesVector, dNodalMassMatrixTrTimesVector

"""
    M = getNodalMassMatrix(mesh::OcTreeMeshFV, sigma::Vector)

Compute finite volume nodal mass matrix with coefficient `sigma` on OcTree
mesh `mesh`.

`getNodalMassMatrix` requires a scalar coefficient, that is, `length(sigma)`
 must equal the number of nodes `mesh.nc`.

The first time that one of `getNodalMassMatrix`, `getdNodalMassMatrix`,
`dNodalMassMatrixTimesVector` or `dNodalMassMatrixTrTimesVector` is called
for the mesh, integration weights and data structure are precomputed and
stored internally. Subsequent calls execute much faster by exploiting the
stored data structure.

The mass matrix integration neglects hanging nodes. For best results,
elimiate hanging nodes as follows:

    M = getNodalMassMatrix(mesh, sigma)
    N, = getNodalConstraints(mesh)
    M = N' * M * N
"""
function getNodalMassMatrix(mesh::OcTreeMeshFV,sigma::Vector)
    P = setupNodalMassMatrix(mesh, sigma)
    nzval = P.A * sigma
    M = SparseMatrixCSC(P.n, P.n, P.colptr, P.rowval, nzval)
    return M
end

"""
    dM = getdNodalMassMatrix(mesh::OcTreeMeshFV, sigma::Vector, v::Vector)

Compute directional derivative of nodal mass matrix w.r.t. coefficient `sigma`
```math
\\nabla_\\sigma (M(\\sigma) v)
```

`length(v)` must equal the number of nodes.

`getdNodalMassMatrix` returns of sparse matrix of size
`(length(v), length(sigma))`.

See documentation of `getNodalMassMatrix` for more details.
"""
function getdNodalMassMatrix{T<:Number}(mesh::OcTreeMeshFV, sigma::Vector, v::Vector{T})
    P = setupNodalMassMatrix(mesh, sigma)
    dM = SparseMatrixCSC(P.A.m, P.A.n, copy(P.A.colptr), copy(P.A.rowval), T.(P.A.nzval))
    DiagTimesM(v, dM)
    return dM
end

# for backwards compatibility
getdNodalMassMatrix(mesh::OcTreeMeshFV, v::Vector) = getdNodalMassMatrix(mesh, ones(M.nc), v)

"""
    dMx = dNodalMassMatrixTimesVector(mesh::OcTreeMeshFV, sigma::Vector, v::Vector, x::Vector)

Compute matrix-vector product between directional derivative of nodal mass matrix
 w.r.t. coefficient `sigma`
```math
\\nabla_\\sigma (M(\\sigma) v)
```
and vector `x`.

`dNodalMassMatrixTimesVector(mesh, sigma, v, x)` is a fast implementation of
`getdNodalMassMatrix(mesh, sigma, v) * x`.

`length(v)` must equal the number of nodes;
`x` must satisfy `length(x) == length(sigma)`.

`dNodalMassMatrixTimesVector` returns a vector of size `length(v)`.

See documentation of `getNodalMassMatrix` for more details.
"""
function dNodalMassMatrixTimesVector{T<:Number}(mesh::OcTreeMeshFV, sigma::Vector, v::Vector{T}, x::Vector)
    P = setupNodalMassMatrix(mesh, sigma)
    dMx = T.(P.A * x)
    dMx .*= v
    return dMx
end

"""
    dMTx = dNodalMassMatrixTrTimesVector(mesh::OcTreeMeshFV, sigma::Vector, v::Vector, x::Vector)

Compute matrix-vector product between (complex) transpose of  directional
derivative of nodal mass matrix
 w.r.t. coefficient `sigma`
```math
\\nabla_\\sigma (M(\\sigma) v)
```
and vector `x`.

`dNodalMassMatrixTimesTrVector(mesh, sigma, v, x)` is a fast implementation of
`getdNodalMassMatrix(mesh, sigma, v)' * x`.

`length(v)` must equal the number of nodes;
`x` must satisfy `length(x) == length(v)`.

`dNodalMassMatrixTimesTrVector` returns a vector of size `length(sigma)`.

See documentation of `getNodalMassMatrix` for more details.
"""
function dNodalMassMatrixTrTimesVector{T<:Number}(mesh::OcTreeMeshFV, sigma::Vector, v::Vector{T}, x::Vector{T})
    P = setupNodalMassMatrix(mesh, sigma)
    u = conj.(v) .* x
    dMTx = P.A' * u
    return dMTx
end

### Private (not exported) methods below

"""
    P = setupNodalMassMatrix(mesh, sigma)

Precompute and return stored data structure for nodal mass matrix integration.
"""
function setupNodalMassMatrix(mesh, sigma)

    na = length(sigma)
    nc = mesh.nc
    N  = typeof(nc)

    @assert na == nc "Invalid size of sigma"

    if !haskey(mesh.Pn, na)

        n = mesh.nn
        i, w = getNodalMassMatrixQuadrature(mesh)
        j = repmat(one(N):nc, 8)
        a = repmat(w, 8)
        A = sparse(i, j, a, n, na)
        colptr = collect(N,1:n+1)
        rowval = collect(N,1:n)
        colval = Array{N}(0) # unused

        mesh.Pn[na] = MassMatrix(Int(n), A, rowval, colptr, colval)

    end

    return mesh.Pn[na]

end

"""
    n, w = getNodalMassMatrixQuadrature(mesh)

Compute quadrature points and weights for integration of nodal mass matrix.

`getNodalMassMatrixQuadrature` returns a list of nodes and weights
for eight quadrature points, the corners of each cell.
"""
function getNodalMassMatrixQuadrature(mesh)

# To illustrate the integration, consider the 2D case (QuadTree) and, in
# particular, the cell C1 of size 2 * h which has two right neighbors of
# of size h:
#
#    n3              n4
#      +-------------+------+------+
#      |q3         q4|      |      |
#      |             |      |      |
#      |             |      |      |
#      |      C1     +------+------+
#      |             |      |      |
#      |             |      |      |
#      |q1         q2|      |      |
#      +-------------+------+------+
#    n1              n2
#
# The integration of the mass matrix for the cell is performed in three
# steps:
# 1. Interpolate the nodal function, discretized by n1, n2, n3 and n4,
#    to the four quadrature points qi (i = 1...4)
#      N(qi) = P(n1,n2,n3,n4)
# 2. Compute the inner product weighted by the scalar coefficient C1
#      Ii = N(qi) * C1 * N(qi)   (i = 1...4)
# 3. Integrate using numerical quadrature
#      M1 = (I1 + I2 + I3 + I4) / 4 * (2 * h)^2
#
# We use the following nearest neighbour interpolation
#   N(q1) = n1
#   N(q2) = n2
#   N(q3) = n3
#   N(q4) = n4
#
# The mass matrix has up to 9 non-zero entries per row/column.
#
# Note that this integration method neglects hanging nodes. The
# mass matrix
#   M1 = getNodalMassMatrix(mesh, sigma)
# differs from
#   Ae = getNodalToCellCenterMatrix(mesh) # nodal to cell center average
#   V  = getVolume(mesh)                  # cell volume
#   M2 = sdiag(Ae' * (V * sigma))
# unless the OcTree degenerates to a regular mesh.
#
# C. Schwarzbach, April 2014

# nodal numbering (implies location of nodal)
N = getNodalNumbering(mesh)

# cells
i,j,k,bsz = find3(mesh.S)
mx,my,mz  = size(mesh.S)

# locate nodal numbers for quadrature points
i0 = i
i1 = i + bsz
j0 = j
j1 = j + bsz
k0 = k
k1 = k + bsz

n000 = getNodesFromIndices(N.SV,(mx+1,my+1,mz+1),i0,j0,k0)
n100 = getNodesFromIndices(N.SV,(mx+1,my+1,mz+1),i1,j0,k0)
n010 = getNodesFromIndices(N.SV,(mx+1,my+1,mz+1),i0,j1,k0)
n110 = getNodesFromIndices(N.SV,(mx+1,my+1,mz+1),i1,j1,k0)
n001 = getNodesFromIndices(N.SV,(mx+1,my+1,mz+1),i0,j0,k1)
n101 = getNodesFromIndices(N.SV,(mx+1,my+1,mz+1),i1,j0,k1)
n011 = getNodesFromIndices(N.SV,(mx+1,my+1,mz+1),i0,j1,k1)
n111 = getNodesFromIndices(N.SV,(mx+1,my+1,mz+1),i1,j1,k1)

# collect everything
n = vcat(n000,n100,n010,n110,n001,n101,n011,n111)

# integration weights
w = bsz.^3 * (prod(mesh.h) / 8)

return n, w

end

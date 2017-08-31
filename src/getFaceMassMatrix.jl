export getFaceMassMatrix, getdFaceMassMatrix,
       dFaceMassMatrixTimesVector, dFaceMassMatrixTrTimesVector

"""
    M = getFaceMassMatrix(mesh::OcTreeMeshFV, sigma::Vector)

Compute finite volume face mass matrix with coefficient `sigma` on OcTree
mesh `mesh`.

`getFaceMassMatrix` distinguishes three cases:
 - `length(sigma) == mesh.nc`: scalar coefficient (isotropy)
 - `length(sigma) == mesh.nc * 3`: diagonal tensor coefficient (diagonal anisotropy)
 - `length(sigma) == mesh.nc * 6`: symmetric tensor coefficient (general anisotropy)

The tensor coefficients must be stored in the following order:
 - `sigma = vcat(sigma_xx, sigma_yy, sigma_zz)`
 - `sigma = vcat(sigma_xx, sigma_yy, sigma_zz, sigma_xy, sigma_xz, sigma_yz)`
where the vectors `sigma_xx`, ... contain the respective entries of the
tensor for all OcTree cells.

The first time that one of `getFaceMassMatrix`, `getdFaceMassMatrix`,
`dFaceMassMatrixTimesVector` or `dFaceMassMatrixTrTimesVector` is called for
the mesh and for a particular type of coefficient, integration weights and
data structure are precomputed and stored internally. Subsequent calls for
the same type of coefficient execute much faster by exploiting the stored
data structure.
"""
function getFaceMassMatrix(mesh::OcTreeMeshFV,sigma::Vector)
    P = setupFaceMassMatrix(mesh, sigma)
    nzval = P.A * sigma
    M = SparseMatrixCSC(P.n, P.n, P.colptr, P.rowval, nzval)
    return M
end

"""
    dM = getdFaceMassMatrix(mesh::OcTreeMeshFV, sigma::Vector, v::Vector)

Compute directional derivative of face mass matrix w.r.t. coefficient `sigma`
```math
\\nabla_\\sigma (M(\\sigma) v)
```

`length(v)` must equal the number of faces.

`getdFaceMassMatrix` returns of sparse matrix of size
`(length(v), length(sigma))`.

See documentation of `getFaceMassMatrix` for more details.
"""
function getdFaceMassMatrix{T<:Number}(mesh::OcTreeMeshFV, sigma::Vector, v::Vector{T})
    P = setupFaceMassMatrix(mesh, sigma)
    nz = length(P.colval)
    if nz == 0 # isotropic or diagonally anisotropic
        dM = SparseMatrixCSC(P.A.m, P.A.n, copy(P.A.colptr), copy(P.A.rowval), T.(P.A.nzval))
        DiagTimesM(v, dM)
    else # generally anisotropic
        V = SparseMatrixCSC(P.n, nz, collect(1:nz+1), P.rowval, v[P.colval])
        dM = V * P.A
    end
    return dM
end

# for backwards compatibility
getdFaceMassMatrix(mesh::OcTreeMeshFV, v::Vector) = getdFaceMassMatrix(mesh, ones(M.nc), v)

"""
    dMx = dFaceMassMatrixTimesVector(mesh::OcTreeMeshFV, sigma::Vector, v::Vector, x::Vector)

Compute matrix-vector product between directional derivative of face mass matrix
 w.r.t. coefficient `sigma`
```math
\\nabla_\\sigma (M(\\sigma) v)
```
and vector `x`.

`dFaceMassMatrixTimesVector(mesh, sigma, v, x)` is a fast implementation of
`getdFaceMassMatrix(mesh, sigma, v) * x`.

`length(v)` must equal the number of faces;
`x` must satisfy `length(x) == length(sigma)`.

`dFaceMassMatrixTimesVector` returns a vector of size `length(v)`.

See documentation of `getFaceMassMatrix` for more details.
"""
function dFaceMassMatrixTimesVector{T<:Number}(mesh::OcTreeMeshFV, sigma::Vector, v::Vector{T}, x::Vector)
    P = setupFaceMassMatrix(mesh, sigma)
    nz = length(P.colval)
    if nz == 0 # isotropic or diagonally anisotropic
        dMx = T.(P.A * x)
        dMx .*= v
    else # generally anisotropic
        V = SparseMatrixCSC(P.n, nz, collect(1:nz+1), P.rowval, v[P.colval])
        dMx = V * (P.A * x)
    end
    return dMx
end

"""
    dMTx = dFaceMassMatrixTrTimesVector(mesh::OcTreeMeshFV, sigma::Vector, v::Vector, x::Vector)

Compute matrix-vector product between (complex) transpose of  directional
derivative of face mass matrix
 w.r.t. coefficient `sigma`
```math
\\nabla_\\sigma (M(\\sigma) v)
```
and vector `x`.

`dFaceMassMatrixTimesTrVector(mesh, sigma, v, x)` is a fast implementation of
`getdFaceMassMatrix(mesh, sigma, v)' * x`.

`length(v)` must equal the number of faces;
`x` must satisfy `length(x) == length(v)`.

`dFaceMassMatrixTimesTrVector` returns a vector of size `length(sigma)`.

See documentation of `getFaceMassMatrix` for more details.
"""
function dFaceMassMatrixTrTimesVector{T<:Number}(mesh::OcTreeMeshFV, sigma::Vector, v::Vector{T}, x::Vector{T})
    P = setupFaceMassMatrix(mesh, sigma)
    nz = length(P.colval)
    if nz == 0 # isotropic or diagonally anisotropic
        u = conj.(v) .* x
        dMTx = P.A' * u
    else # generally anisotropic
        V = SparseMatrixCSC(P.n, nz, collect(1:nz+1), P.rowval, v[P.colval])
        u = V' * x
        dMTx = P.A' * u
    end
    return dMTx
end

### Private (not exported) methods below

"""
    P = setupFaceMassMatrix(mesh, sigma)

Precompute and return stored data structure for face mass matrix integration
for isotropic, diagonally anisotropic or generally anisotropic coefficient.
"""
function setupFaceMassMatrix(mesh, sigma)

    na = length(sigma)
    nc = mesh.nc

    @assert in(na, [nc, 3*nc, 6*nc]) "Invalid size of sigma"

    if !haskey(mesh.Pf, na)

        n = sum(mesh.nf)
        fx, fy, fz, w = getFaceMassMatrixQuadrature(mesh)

        if na == nc

            i = vcat(fx, fy, fz)
            j = repmat(1:nc, 24)
            a = repmat(w, 24)
            A = sparse(i, j, a, n, na)
            colptr = collect(1:n+1)
            rowval = collect(1:n)
            colval = Array{Int64}(0) # unused

        elseif na == 3 * nc

            i = vcat(fx, fy, fz)
            j = vcat(repmat(     1:  nc, 8),
                     repmat(  nc+1:2*nc, 8),
                     repmat(2*nc+1:3*nc, 8))
            a = repmat(w, 24)
            A = sparse(i, j, a, n, na)
            colptr = collect(1:n+1)
            rowval = collect(1:n)
            colval = Array{Int64}(0) # unused

        elseif na == 6 * nc

            # find nonzero pattern of mass matrix
            i = vcat(fx, fy, fz, fx, fx, fy, fy, fz, fz)
            j = vcat(fx, fy, fz, fy, fz, fz, fx, fx, fy)
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
            j = vcat(repmat(     1:  nc, 8),
                     repmat(  nc+1:2*nc, 8),
                     repmat(2*nc+1:3*nc, 8),
                     repmat(3*nc+1:4*nc, 8),
                     repmat(4*nc+1:5*nc, 8),
                     repmat(5*nc+1:6*nc, 8),
                     repmat(3*nc+1:4*nc, 8),
                     repmat(4*nc+1:5*nc, 8),
                     repmat(5*nc+1:6*nc, 8))
            a = repmat(w, 72)
            A = sparse(i, j, a, nz, na)

        end

        mesh.Pf[na] = MassMatrix(n, A, rowval, colptr, colval)

    end

    return mesh.Pf[na]

end

"""
    fx, fy, fz, w = getFaceMassMatrixQuadrature(mesh)

Compute quadrature points and weights for integration of face mass matrix.

`getFaceMassMatrixQuadrature` returns a list of faces and weights
for eight quadrature points, the corners of each cell.
"""
function getFaceMassMatrixQuadrature(mesh)

# To illustrate the integration, consider the 2D case (QuadTree) and, in
# particular, the cell C1 of size 2 * h which has two right neighbors of
# of size h:
#
#            fx2
#      +-------------+------+------+
#      |q3         q4|      |      |
#      |             | fy3  |      |
#      |             |      |      |
#  fy1 |      C1     +------+------+
#      |             |      |      |
#      |             | fy2  |      |
#      |q1         q2|      |      |
#      +-------------+------+------+
#            fx1
#
# The integration of the mass matrix for the cell is performed in three
# steps:
# 1. Interpolate the face function, which contains the x- and y-field
#    components fx1, fx2, fy1, fy2, fy3 staggered on the x- and y- faces,
#    to the four quadrature points qi (i = 1...4)
#      Fx(qi) = Px(fx1,fx2)
#      Fy(qi) = Py(fy1,fy2,fy3)
# 2. Compute the inner product of the collocated field, weighted by the
#    tensor C1
#      I1 = (Fx(qi), Fy(qi)) * C1 * (Fx(qi), Fy(qi))^T   (i = 1...4)
# 3. Integrate using numerical quadrature
#      M1 = (I1 + I2 + I3 + I4) / 4 * (2 * h)^2
#
# We use the following nearest neighbour interpolation
#   Fx(q1) = fx1, Fy(q1) = fy1
#   Fx(q2) = fx1, Fy(q2) = fy2
#   Fx(q3) = fx2, Fy(q3) = fy1
#   Fx(q4) = fx2, Fy(q4) = fy3
#
# The mass matrix has up to 17 non-zero entries per row/column.
#
# C. Schwarzbach, April 2014

# face numbering (implies location of face)
FX,FY,FZ  = getFaceNumbering(mesh)
nfx = nnz(FX)
nfy = nnz(FY)
nfz = nnz(FZ)

# cells
i,j,k,bsz = find3(mesh.S)
mx,my,mz  = size(mesh.S)

# locate face numbers for quadrature points
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

# x-faces
fx000 = getNodesFromIndices(FX.SV,(mx+1,my,mz),i0,j0,k0)
fx100 = getNodesFromIndices(FX.SV,(mx+1,my,mz),i1,j0,k0)
fx010 = getNodesFromIndices(FX.SV,(mx+1,my,mz),i0,j2,k0)
fx110 = getNodesFromIndices(FX.SV,(mx+1,my,mz),i1,j2,k0)
fx001 = getNodesFromIndices(FX.SV,(mx+1,my,mz),i0,j0,k2)
fx101 = getNodesFromIndices(FX.SV,(mx+1,my,mz),i1,j0,k2)
fx011 = getNodesFromIndices(FX.SV,(mx+1,my,mz),i0,j2,k2)
fx111 = getNodesFromIndices(FX.SV,(mx+1,my,mz),i1,j2,k2)

merge!(fx010, fx000)
merge!(fx001, fx000)
merge!(fx011, fx000)
merge!(fx110, fx100)
merge!(fx101, fx100)
merge!(fx111, fx100)

# y-faces
fy000 = getNodesFromIndices(FY.SV,(mx,my+1,mz),i0,j0,k0)
fy100 = getNodesFromIndices(FY.SV,(mx,my+1,mz),i2,j0,k0)
fy010 = getNodesFromIndices(FY.SV,(mx,my+1,mz),i0,j1,k0)
fy110 = getNodesFromIndices(FY.SV,(mx,my+1,mz),i2,j1,k0)
fy001 = getNodesFromIndices(FY.SV,(mx,my+1,mz),i0,j0,k2)
fy101 = getNodesFromIndices(FY.SV,(mx,my+1,mz),i2,j0,k2)
fy011 = getNodesFromIndices(FY.SV,(mx,my+1,mz),i0,j1,k2)
fy111 = getNodesFromIndices(FY.SV,(mx,my+1,mz),i2,j1,k2)

merge!(fy100, fy000)
merge!(fy001, fy000)
merge!(fy101, fy000)
merge!(fy110, fy010)
merge!(fy011, fy010)
merge!(fy111, fy010)

# z-faces
fz000 = getNodesFromIndices(FZ.SV,(mx,my,mz+1),i0,j0,k0)
fz100 = getNodesFromIndices(FZ.SV,(mx,my,mz+1),i2,j0,k0)
fz010 = getNodesFromIndices(FZ.SV,(mx,my,mz+1),i0,j2,k0)
fz110 = getNodesFromIndices(FZ.SV,(mx,my,mz+1),i2,j2,k0)
fz001 = getNodesFromIndices(FZ.SV,(mx,my,mz+1),i0,j0,k1)
fz101 = getNodesFromIndices(FZ.SV,(mx,my,mz+1),i2,j0,k1)
fz011 = getNodesFromIndices(FZ.SV,(mx,my,mz+1),i0,j2,k1)
fz111 = getNodesFromIndices(FZ.SV,(mx,my,mz+1),i2,j2,k1)

merge!(fz100, fz000)
merge!(fz010, fz000)
merge!(fz110, fz000)
merge!(fz101, fz001)
merge!(fz011, fz001)
merge!(fz111, fz001)

# collect everything
fx = vcat(fx000,fx100,fx010,fx110,fx001,fx101,fx011,fx111)
fy = vcat(fy000,fy100,fy010,fy110,fy001,fy101,fy011,fy111)
fz = vcat(fz000,fz100,fz010,fz110,fz001,fz101,fz011,fz111)

# shift fy and fz
fy .+= nfx
fz .+= (nfx + nfy)

# integration weights
w = bsz.^3 * (prod(mesh.h) / 8)

return fx, fy, fz, w

end

for Tn2 in intTypes
for Tn in intTypes
@testset "Constraint matrices with intType1=$Tn, intType2=$Tn2" begin
using JOcTree
using Base.Test

include("getFunction.jl")

# Test the following matrices:
#   getNodalGradientMatrix, getCurlMatrix, getDivergenceMatrix
# These routines (mostly) depend on:
#   getFaceSizeNumbering, getEdgeSizeNumbering, getCellNumbering,
#   getNodalNumbering, getNumberOfNeighbors, sparse3, etc.


println("Testing: getNodalGradientMatrix, getCurlMatrix, getDivergenceMatrix")

#-------------------------------------------

n = (64, 64, 32)   # underlying mesh
h = 1 ./ collect(n)   # cell size
x0 = [0. , 0. , 0.]

nrand = 5
S = randomOctreeMesh(Tn, Tn2, n, nrand )

M = getOcTreeMeshFV(S, h, x0=x0)

# exportUBCOcTreeMesh("mesh.txt", M)

#-------------------------------------------

Ne,Qe, = getEdgeConstraints(M)
Nf,Qf, = getFaceConstraints(M)
Nn,Qn, = getNodalConstraints(M)

GRAD = getNodalGradientMatrix(M)
GRAD = Qe*GRAD*Nn
CURL = getCurlMatrix(M)
CURL = Qf*CURL*Ne
DIV = getDivergenceMatrix(M)
DIV = DIV*Nf

DIVCURL = DIV * CURL
dropzeros!(DIVCURL)
CURLGRAD = CURL * GRAD
dropzeros!(CURLGRAD)

println()
println("nnz(DIVCURL)   ", nnz(DIVCURL))
println("nnz(CURLGRAD)  ", nnz(CURLGRAD))

@test nnz(DIVCURL) == 0
@test nnz(CURLGRAD) == 0

println()

# Check that constraint matrices don't do anything for a uniform mesh
n = 16
SV = sparsevec(1:n^3,1,n^3)
S  = SparseArray3D(SV,(n,n,n))
h  = fill(1/n,3)
M2 = getOcTreeMeshFV(S,h)

Ne,Qe, = getEdgeConstraints(M2)
Nf,Qf, = getFaceConstraints(M2)
Nn,Qn, = getNodalConstraints(M2)

GRAD  = getNodalGradientMatrix(M2)
GRADc = Qe*GRAD*Nn
CURL  = getCurlMatrix(M2)
CURLc = Qf*CURL*Ne
DIV   = getDivergenceMatrix(M2)
DIVc  = DIV*Nf
@test GRAD == GRADc
@test CURL == CURLc
@test DIV == DIVc
end
end
end

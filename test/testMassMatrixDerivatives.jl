println("Testing: getdEdgeMassMatrix, getdFaceMassMatrix, dEdgeMassMatrixTimesVector, dEdgeMassMatrixTrTimesVector")

@testset "Mass matrix derivatives" begin

using jInv.Utils

# Get random mesh
n = [128 , 128 , 128]   # underlying mesh
h = 1 ./ n   # cell size
x0 = [0. , 0. , 0.]

nrand = 5
S = randomOctreeMesh( n, nrand )
M = getOcTreeMeshFV(S, h, x0=x0)

# Get random edge and face fields
ne     = sum(M.ne)
efield = randn(ne)
nf     = sum(M.nf)
ffield = randn(nf)

# Get random cell-centered variables, one for each
# level of anisotropy
nc   = M.nc
sig1 = rand(nc)
sig2 = rand(3*nc)
sig3 = rand(6*nc)

# Edges
Me(x::Vector)  = getEdgeMassMatrix(M,x) * efield
dMe(x::Vector) = getdEdgeMassMatrix(M,x,efield)

p1, = checkDerivative(Me,dMe,sig1)
p2, = checkDerivative(Me,dMe,sig2)
p3, = checkDerivative(Me,dMe,sig3)
@test p1
@test p2
@test p3

# Faces
Mf(x::Vector)  = getFaceMassMatrix(M,x) * ffield
dMf(x::Vector) = getdFaceMassMatrix(M,x,ffield)

p4, = checkDerivative(Mf,dMf,sig1)
p5, = checkDerivative(Mf,dMf,sig2)
p6, = checkDerivative(Mf,dMf,sig3)

@test p4
@test p5
@test p6

# Derivative of edge mass matrix times vector products
dsig1 = rand(nc)
dsig2 = rand(3*nc)
dsig3 = rand(6*nc)
defield = rand(ne)

dMe1 = getdEdgeMassMatrix(M,sig1,efield)
dMe2 = getdEdgeMassMatrix(M,sig2,efield)
dMe3 = getdEdgeMassMatrix(M,sig3,efield)

dMe1dsig1de = dot(defield, dMe1 * dsig1)
dMe2dsig2de = dot(defield, dMe2 * dsig2)
dMe3dsig3de = dot(defield, dMe3 * dsig3)

dMe1dsig1 = dEdgeMassMatrixTimesVector(M,sig1,efield,dsig1)
dMe2dsig2 = dEdgeMassMatrixTimesVector(M,sig2,efield,dsig2)
dMe3dsig3 = dEdgeMassMatrixTimesVector(M,sig3,efield,dsig3)

dMe1de = dEdgeMassMatrixTrTimesVector(M,sig1,efield,defield)
dMe2de = dEdgeMassMatrixTrTimesVector(M,sig2,efield,defield)
dMe3de = dEdgeMassMatrixTrTimesVector(M,sig3,efield,defield)

@test dot(defield, dMe1dsig1) ≈ dMe1dsig1de
@test dot(defield, dMe2dsig2) ≈ dMe2dsig2de
@test dot(defield, dMe3dsig3) ≈ dMe3dsig3de

@test dot(dsig1, dMe1de) ≈ dMe1dsig1de
@test dot(dsig2, dMe2de) ≈ dMe2dsig2de
@test dot(dsig3, dMe3de) ≈ dMe3dsig3de

@test dot(defield, dMe1dsig1) ≈ dot(dsig1, dMe1de)
@test dot(defield, dMe2dsig2) ≈ dot(dsig2, dMe2de)
@test dot(defield, dMe3dsig3) ≈ dot(dsig3, dMe3de)

end
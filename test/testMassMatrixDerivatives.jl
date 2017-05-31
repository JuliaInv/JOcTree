println("Testing: getdEdgeMassMatrix, getdFaceMassMatrix")

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
Me1(x::Vector) = getEdgeMassMatrix(M,x)*efield
Me2(x::Vector) = getEdgeMassMatrix(M,x)*efield
Me3(x::Vector) = getEdgeMassMatrix(M,x)*efield


dMe1(x::Vector) = getdEdgeMassMatrix(M,sig1,efield)
dMe2(x::Vector) = getdEdgeMassMatrix(M,sig2,efield)
dMe3(x::Vector) = getdEdgeMassMatrix(M,sig3,efield)

p1, = checkDerivative(Me1,dMe1,sig1)
p2, = checkDerivative(Me2,dMe2,sig2)
p3, = checkDerivative(Me3,dMe3,sig3)
@test p1
@test p2
@test p3

#-------------------------------------------

# Faces
Mf1(x::Vector) = getFaceMassMatrix(M,x)*ffield
Mf2(x::Vector) = getFaceMassMatrix(M,x)*ffield
Mf3(x::Vector) = getFaceMassMatrix(M,x)*ffield

dMf1(x::Vector) = getdFaceMassMatrix(M,sig1,ffield)
dMf2(x::Vector) = getdFaceMassMatrix(M,sig2,ffield)
dMf3(x::Vector) = getdFaceMassMatrix(M,sig3,ffield)

p4, = checkDerivative(Mf1,dMf1,sig1)
p5, = checkDerivative(Mf2,dMf2,sig2)
p6, = checkDerivative(Mf3,dMf3,sig3)

@test p4
@test p5
@test p6
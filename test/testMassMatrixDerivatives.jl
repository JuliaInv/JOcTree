println("Testing: getdNodalMassMatrix, dNodalMassMatrixTimesVector, dNodalMassMatrixTrTimesVector,")
println("         getdEdgeMassMatrix, dEdgeMassMatrixTimesVector, dEdgeMassMatrixTrTimesVector,")
println("         getdFaceMassMatrix, dFaceMassMatrixTimesVector, dFaceMassMatrixTrTimesVector")

@testset "Mass matrix derivatives" begin

using jInv.Utils

# Get random mesh
n = [128 , 128 , 128]   # underlying mesh
h = 1 ./ n   # cell size
x0 = [0. , 0. , 0.]

nrand = 5
S = randomOctreeMesh( n, nrand )
mesh = getOcTreeMeshFV(S, h, x0=x0)

# Get random nodal, edge and face fields
nn       = mesh.nn
nreal    = randn(nn)
ncomplex = complex.(nreal, randn(nn))
ne       = sum(mesh.ne)
ereal    = randn(ne)
ecomplex = complex.(ereal, randn(ne))
nf       = sum(mesh.nf)
freal    = randn(nf)
fcomplex = complex.(freal, randn(nf))

# Get random cell-centered variables, one for each
# level of anisotropy
nc   = mesh.nc
sig1 = rand(nc)
sig2 = rand(3*nc)
sig3 = rand(6*nc)

# Nodes
for u in (nreal, ncomplex)
  M(x::Vector)  = getNodalMassMatrix(mesh, x) * u
  dM(x::Vector) = getdNodalMassMatrix(mesh, x, u)
  for sig in (sig1,)
    p, = checkDerivative(M, dM, sig)
    @test p
  end
end

# Edges
for u in (ereal, ecomplex)
  M(x::Vector)  = getEdgeMassMatrix(mesh, x) * u
  dM(x::Vector) = getdEdgeMassMatrix(mesh, x, u)
  for sig in (sig1, sig2, sig3)
    p, = checkDerivative(M, dM, sig)
    @test p
  end
end

# Faces
for u in (freal, fcomplex)
  M(x::Vector)  = getFaceMassMatrix(mesh, x) * u
  dM(x::Vector) = getdFaceMassMatrix(mesh, x, u)
  for sig in (sig1, sig2, sig3)
    p, = checkDerivative(M, dM, sig)
    @test p
  end
end

# Derivative of mass matrix times vector products
dsig1 = rand(nc)
dsig2 = rand(3*nc)
dsig3 = rand(6*nc)
dnreal = rand(nn)
dncomplex = complex.(dnreal, rand(nn))
dereal = rand(ne)
decomplex = complex.(dereal, rand(ne))
dfreal = rand(nf)
dfcomplex = complex.(dfreal, rand(nf))

# Nodes
for (u, du) in zip((nreal, ncomplex), (dnreal, dncomplex))
  for (sig, dsig) in zip((sig1,), (dsig1,))
    dM = getdNodalMassMatrix(mesh, sig, u)
    dMdsigdu = dot(dM * dsig, du)
    dMdsig = dNodalMassMatrixTimesVector(mesh, sig, u, dsig)
    dMdu = dNodalMassMatrixTrTimesVector(mesh, sig, u, du)
    @test dot(dMdsig, du) ≈ dMdsigdu
    @test dot(dsig, dMdu) ≈ dMdsigdu
    @test dot(dMdsig, du) ≈ dot(dsig, dMdu)
  end
end

# Edges
for (u, du) in zip((ereal, ecomplex), (dereal, decomplex))
  for (sig, dsig) in zip((sig1, sig2, sig3), (dsig1, dsig2, dsig3))
    dM = getdEdgeMassMatrix(mesh, sig, u)
    dMdsigdu = dot(dM * dsig, du)
    dMdsig = dEdgeMassMatrixTimesVector(mesh, sig, u, dsig)
    dMdu = dEdgeMassMatrixTrTimesVector(mesh, sig, u, du)
    @test dot(dMdsig, du) ≈ dMdsigdu
    @test dot(dsig, dMdu) ≈ dMdsigdu
    @test dot(dMdsig, du) ≈ dot(dsig, dMdu)
  end
end

# Faces
for (u, du) in zip((freal, fcomplex), (dfreal, dfcomplex))
  for (sig, dsig) in zip((sig1, sig2, sig3), (dsig1, dsig2, dsig3))
    dM = getdFaceMassMatrix(mesh, sig, u)
    dMdsigdu = dot(dM * dsig, du)
    dMdsig = dFaceMassMatrixTimesVector(mesh, sig, u, dsig)
    dMdu = dFaceMassMatrixTrTimesVector(mesh, sig, u, du)
    @test dot(dMdsig, du) ≈ dMdsigdu
    @test dot(dsig, dMdu) ≈ dMdsigdu
    @test dot(dMdsig, du) ≈ dot(dsig, dMdu)
  end
end

end # end of testset
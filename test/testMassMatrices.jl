println("Testing: getNodalMassMatrix, getEdgeMassMatrix, getFaceMassMatrix")

@testset "Mass matrices" begin

# Get random mesh
n = [256, 128, 64]
L = exp.(randn(3))
h = L ./ n
x0 = randn(3)
xn = x0 + L
nrand = 5 # local refinement at nrand random locations
S = randomOctreeMesh(n, nrand)
mesh = getOcTreeMeshFV(S, h, x0=x0)

### Check exact integration of constant and linear coefficient
un = ones(mesh.nn)
uex = vcat(ones(mesh.ne[1]), zeros(mesh.ne[2]), zeros(mesh.ne[3]))
uey = vcat(zeros(mesh.ne[1]), ones(mesh.ne[2]), zeros(mesh.ne[3]))
uez = vcat(zeros(mesh.ne[1]), zeros(mesh.ne[2]), ones(mesh.ne[3]))
ufx = vcat(ones(mesh.nf[1]), zeros(mesh.nf[2]), zeros(mesh.nf[3]))
ufy = vcat(zeros(mesh.nf[1]), ones(mesh.nf[2]), zeros(mesh.nf[3]))
ufz = vcat(zeros(mesh.nf[1]), zeros(mesh.nf[2]), ones(mesh.nf[3]))

V = sum(getVolumeVector(mesh))
a = V * 8 / prod(xn .* xn - x0 .* x0) # normalize c = a * x * y * z such that integral equals volume
cc = ones(mesh.nc) # constant: c = 1
cl = vec(prod(getCellCenteredGrid(mesh), 2)) .* a # linear: c = a * x * y * z

for c1 in (cc, cl)

  c3 = repmat(c1, 3)
  c6 = repmat(c3, 2)

  # nodes
  Mn = getNodalMassMatrix(mesh, c1)
  @test dot(un, Mn * un) ≈ V

  # edges
  Me1 = getEdgeMassMatrix(mesh, c1)
  @test dot(uex, Me1 * uex) ≈ V
  @test dot(uey, Me1 * uey) ≈ V
  @test dot(uez, Me1 * uez) ≈ V

  Me3 = getEdgeMassMatrix(mesh, c3)
  @test dot(uex, Me3 * uex) ≈ V
  @test dot(uey, Me3 * uey) ≈ V
  @test dot(uez, Me3 * uez) ≈ V

  Me6 = getEdgeMassMatrix(mesh, c6)
  @test dot(uex, Me6 * uex) ≈ V
  @test dot(uex, Me6 * uey) ≈ V
  @test dot(uex, Me6 * uez) ≈ V
  @test dot(uey, Me6 * uex) ≈ V
  @test dot(uey, Me6 * uey) ≈ V
  @test dot(uey, Me6 * uez) ≈ V
  @test dot(uez, Me6 * uex) ≈ V
  @test dot(uez, Me6 * uey) ≈ V
  @test dot(uez, Me6 * uez) ≈ V

  # faces
  Mf1 = getFaceMassMatrix(mesh, c1)
  @test dot(ufx, Mf1 * ufx) ≈ V
  @test dot(ufy, Mf1 * ufy) ≈ V
  @test dot(ufz, Mf1 * ufz) ≈ V

  Mf3 = getFaceMassMatrix(mesh, c3)
  @test dot(ufx, Mf3 * ufx) ≈ V
  @test dot(ufy, Mf3 * ufy) ≈ V
  @test dot(ufz, Mf3 * ufz) ≈ V

  Mf6 = getFaceMassMatrix(mesh, c6)
  @test dot(ufx, Mf6 * ufx) ≈ V
  @test dot(ufx, Mf6 * ufy) ≈ V
  @test dot(ufx, Mf6 * ufz) ≈ V
  @test dot(ufy, Mf6 * ufx) ≈ V
  @test dot(ufy, Mf6 * ufy) ≈ V
  @test dot(ufy, Mf6 * ufz) ≈ V
  @test dot(ufz, Mf6 * ufx) ≈ V
  @test dot(ufz, Mf6 * ufy) ≈ V
  @test dot(ufz, Mf6 * ufz) ≈ V

end

# hanging nodes, edges, faces elimination
Nn,Qn, = getNodalConstraints(mesh)
Ne,Qe, = getEdgeConstraints(mesh)
Nf,Qf, = getFaceConstraints(mesh)

un = Qn * un
uex = Qe * uex
uey = Qe * uey
uez = Qe * uez
ufx = Qf * ufx
ufy = Qf * ufy
ufz = Qf * ufz

for c1 in (cc, cl)

  c3 = repmat(c1, 3)
  c6 = repmat(c3, 2)

  # nodes
  Mn = getNodalMassMatrix(mesh, c1)
  Mn = Nn' * Mn * Nn
  @test dot(un, Mn * un) ≈ V

  # edges
  Me1 = getEdgeMassMatrix(mesh, c1)
  Me1 = Ne' * Me1 * Ne
  @test dot(uex, Me1 * uex) ≈ V
  @test dot(uey, Me1 * uey) ≈ V
  @test dot(uez, Me1 * uez) ≈ V

  Me3 = getEdgeMassMatrix(mesh, c3)
  Me3 = Ne' * Me3 * Ne
  @test dot(uex, Me3 * uex) ≈ V
  @test dot(uey, Me3 * uey) ≈ V
  @test dot(uez, Me3 * uez) ≈ V

  Me6 = getEdgeMassMatrix(mesh, c6)
  Me6 = Ne' * Me6 * Ne
  @test dot(uex, Me6 * uex) ≈ V
  @test dot(uex, Me6 * uey) ≈ V
  @test dot(uex, Me6 * uez) ≈ V
  @test dot(uey, Me6 * uex) ≈ V
  @test dot(uey, Me6 * uey) ≈ V
  @test dot(uey, Me6 * uez) ≈ V
  @test dot(uez, Me6 * uex) ≈ V
  @test dot(uez, Me6 * uey) ≈ V
  @test dot(uez, Me6 * uez) ≈ V

  # faces
  Mf1 = getFaceMassMatrix(mesh, c1)
  Mf1 = Nf' * Mf1 * Nf
  @test dot(ufx, Mf1 * ufx) ≈ V
  @test dot(ufy, Mf1 * ufy) ≈ V
  @test dot(ufz, Mf1 * ufz) ≈ V

  Mf3 = getFaceMassMatrix(mesh, c3)
  Mf3 = Nf' * Mf3 * Nf
  @test dot(ufx, Mf3 * ufx) ≈ V
  @test dot(ufy, Mf3 * ufy) ≈ V
  @test dot(ufz, Mf3 * ufz) ≈ V

  Mf6 = getFaceMassMatrix(mesh, c6)
  Mf6 = Nf' * Mf6 * Nf
  @test dot(ufx, Mf6 * ufx) ≈ V
  @test dot(ufx, Mf6 * ufy) ≈ V
  @test dot(ufx, Mf6 * ufz) ≈ V
  @test dot(ufy, Mf6 * ufx) ≈ V
  @test dot(ufy, Mf6 * ufy) ≈ V
  @test dot(ufy, Mf6 * ufz) ≈ V
  @test dot(ufz, Mf6 * ufx) ≈ V
  @test dot(ufz, Mf6 * ufy) ≈ V
  @test dot(ufz, Mf6 * ufz) ≈ V

end

### Consistency checks for anisotropy

## isotropic
c1 = rand(mesh.nc)
c3 = repmat(c1, 3)
c6 = vcat(c3, zeros(3 * mesh.nc))

# edges
Me1 = getEdgeMassMatrix(mesh, c1)
Me3 = getEdgeMassMatrix(mesh, c3)
Me6 = getEdgeMassMatrix(mesh, c6)
@test Me1 ≈ Me3
@test Me1 ≈ Me6

# faces
Mf1 = getFaceMassMatrix(mesh, c1)
Mf3 = getFaceMassMatrix(mesh, c3)
Mf6 = getFaceMassMatrix(mesh, c6)
@test Mf1 ≈ Mf3
@test Mf1 ≈ Mf6

## diagonally anisotropic
c3 = rand(3 * mesh.nc)
c6 = vcat(c3, zeros(3 * mesh.nc))

# edges
Me3 = getEdgeMassMatrix(mesh, c3)
Me6 = getEdgeMassMatrix(mesh, c6)
@test Me3 ≈ Me6

# faces
Mf3 = getFaceMassMatrix(mesh, c3)
Mf6 = getFaceMassMatrix(mesh, c6)
@test Mf3 ≈ Mf6

end # end of testset
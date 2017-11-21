println("Testing: getInterpolationMatrix")
@testset "Mesh interpolation" begin
for Tn2 in intTypes
for Tn in intTypes
@testset "Mesh interpolation with intType1=$Tn, intType2=$Tn2" begin

using jInv.Utils

# Get random mesh
n = (128 , 128 , 128)   # underlying mesh
h = 1 ./ collect(n)   # cell size
x0 = [0. , 0. , 0.]

nrand = 5
S1 = randomOctreeMesh(Tn, Tn2, n, nrand )
M1 = getOcTreeMeshFV(S1, h, x0=x0)

S2 = createOcTreeFromBox(x0[1],x0[2],x0[3],n[1],n[2],n[3],h[1],h[2],h[3],
                         0.5,0.5,0.5,0.5,0.5,0.5,2,2,1,Tn,Tn2)
M2 = getOcTreeMeshFV(S2, h, x0=x0)

P12 = getInterpolationMatrix(M1,M2)
P21 = getInterpolationMatrix(M2,M1)

u1 = ones(M1.nc)
u2 = ones(M2.nc)

@test all(P12 * u1 .== 1.0)
@test all(P21 * u2 .== 1.0)

# exact mapping if N >= largest block size in S1 and S2
N = round(Int, log2(max(maximum(nonzeros(S1)), maximum(nonzeros(S2)))))
PN12 = getInterpolationMatrix(M1,M2,N=N)
PN21 = getInterpolationMatrix(M2,M1,N=N)

@test P12 == PN12
@test P21 == PN21

end
end
end
end

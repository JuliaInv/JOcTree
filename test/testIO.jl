println("Testing: exportUBCOcTreeMesh, importUBCOcTreeMesh, exportUBCOcTreeModel, importUBCOcTreeModel")

@testset "IO" begin

n = [128, 128, 128]
h = rand(1.:100.0,3)   # cell size
x0 = rand(1.:100.0,3)

nrand = 5
S = randomOctreeMesh(n, nrand)

meshOut = getOcTreeMeshFV(S, h; x0=x0)
exportUBCOcTreeMesh("mesh.msh", meshOut)
meshIn = importUBCOcTreeMesh("mesh.msh")
@test meshIn==meshOut

for modelOut in [
    rand(1:0.1:100, meshOut.nc),
    rand(1:0.1:100, meshOut.nc,3),
    rand(1:0.1:100, meshOut.nc,6)]
    
    exportUBCOcTreeModel("model.mod", meshOut, modelOut)
    modelIn = importUBCOcTreeModel("model.mod", meshIn)
    @test modelIn==modelOut  
    
end

for modelOut in [
    rand(Bool, meshOut.nc),
    rand(Bool, meshOut.nc, 3),
    rand(1:100, meshOut.nc),
    rand(1:100, meshOut.nc, 3),
    rand(1:0.1:100, meshOut.nc),
    rand(1:0.1:100, meshOut.nc, 3)]
    
    exportUBCOcTreeModel("model.mod", meshOut, modelOut)
    modelIn = importUBCOcTreeModel("model.mod", meshIn, eltype(modelOut))
    @test modelIn==modelOut
    @test eltype(modelIn)==eltype(modelOut)
    
end

rm("mesh.msh")
rm("model.mod")

end

include("randomOctreeMesh.jl") 
S = randomOctreeMesh( [256, 256, 256], 5 )


meshOut = getOcTreeMeshFV(S, rand(1.:100.0,3); x0=rand(1.:100.0,3))
modelOut = rand(1:0.1:100, meshOut.nc)
exportUBCOcTreeMesh("mesh.msh", meshOut)
exportUBCOcTreeModel("model.mod", meshOut, modelOut)

meshIn = importUBCOcTreeMesh("mesh.msh")
modelIn = importUBCOcTreeModel("model.mod", meshIn);

@test meshIn==meshOut
@test modelIn==modelOut

rm("mesh.msh")
rm("model.mod")

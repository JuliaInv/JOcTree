
S = initializeOctree([256, 256, 256])
for i in 1:4
    refineInd = unique(rand(1:nnz(S),5))
    S = splitCells(S, refineInd)
end
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

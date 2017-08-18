import JOcTree.hasHDF5

println("Testing: exportUBCOcTreeMesh, importUBCOcTreeMesh, exportUBCOcTreeModel, importUBCOcTreeModel")
if hasHDF5
  println("Testing: exportHDF5OcTreeMesh, importHDF5OcTreeMesh")
end

@testset "IO" begin

n = [128, 128, 128]
h = rand(1.:100.0,3)   # cell size
x0 = rand(1.:100.0,3)

nrand = 5
S  = randomOctreeMesh(n, nrand)
S2 = randomOctreeMesh(n, nrand)
S3 = randomOctreeMesh(n, nrand)

meshOut = getOcTreeMeshFV(S, h; x0=x0)

exportUBCOcTreeMesh("mesh.msh", meshOut)
meshInUBC = importUBCOcTreeMesh("mesh.msh")
@test meshInUBC==meshOut

if hasHDF5
  
  exportHDF5OcTreeMesh("meshHDF5.msh",meshOut)
  meshInHDF5 = importHDF5OcTreeMesh("meshHDF5.msh")
  @test meshInHDF5==meshOut

  # Test storing multiple meshes in one HDF5 file
  m1 = getOcTreeMeshFV( S, h; x0=x0)
  m2 = getOcTreeMeshFV(S2, h; x0=x0)
  m3 = getOcTreeMeshFV(S3, h; x0=x0)
  meshVec = [m1;m2;m3]
  meshDict = Dict([("m1",m1);("m2",m2);("m3",m3)])
  
  exportHDF5OcTreeMesh("vecinHDF5.msh",meshVec)
  exportHDF5OcTreeMesh("dictinHDF5.msh",meshDict)
  
  vecIn = importHDF5OcTreeMesh("vecinHDF5.msh")
  dictIn = importHDF5OcTreeMesh("dictinHDF5.msh")
  for (id,mesh) in zip(["1";"2";"3"],meshVec)
    @test vecIn[id] == mesh
  end
  for (id,mesh) in meshDict
    @test dictIn[id] == mesh
  end

  vecIn = importHDF5OcTreeMesh("vecinHDF5.msh",[1,2,3])
  dictIn = importHDF5OcTreeMesh("dictinHDF5.msh",["m1","m2","m3"])
  for (id,mesh) in zip(["1";"2";"3"],meshVec)
    @test vecIn[id] == mesh
  end
  for (id,mesh) in meshDict
    @test dictIn[id] == mesh
  end
  
  for id = 1:3
    meshInHDF5 = importHDF5OcTreeMesh("vecinHDF5.msh",id)
    @test meshInHDF5 == meshVec[id]
  end
  for id in ("m1","m2","m3")
    meshInHDF5 = importHDF5OcTreeMesh("dictinHDF5.msh",id)
    @test meshInHDF5 == meshDict[id]
  end
  
  @test_throws ErrorException importHDF5OcTreeMesh("dictinHDF5.msh",1)  
  @test_throws ErrorException importHDF5OcTreeMesh("dictinHDF5.msh",[1,2,3])  

end

for modelOut in [
    rand(1:0.1:100, meshOut.nc),
    rand(1:0.1:100, meshOut.nc,3),
    rand(1:0.1:100, meshOut.nc,6)]

    exportUBCOcTreeModel("model.mod", meshOut, modelOut)
    modelIn = importUBCOcTreeModel("model.mod", meshInUBC)
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
    modelIn = importUBCOcTreeModel("model.mod", meshInUBC, eltype(modelOut))
    @test modelIn==modelOut
    @test eltype(modelIn)==eltype(modelOut)

end

rm("mesh.msh")
rm("model.mod")
if hasHDF5
  rm("vecinHDF5.msh")
  rm("dictinHDF5.msh")
  rm("meshHDF5.msh")
end

end

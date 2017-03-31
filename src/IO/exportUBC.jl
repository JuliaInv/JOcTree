
export exportUBCOcTreeMesh, exportUBCOcTreeModel, outputOctreeMesh


"""
    exportUBCOcTreeMesh(fname, mesh)
    
    Writes an OcTree mesh to disk in UBC format.

    Input:
    
        name::AbstractString - File to write
        mesh::OcTreeMesh     - The mesh to export
            
"""
function exportUBCOcTreeMesh(fname::AbstractString, mesh::OcTreeMesh)
# Export OcTree for use with Roman's Fortran codes.

   m1,m2,m3     = mesh.n
   i1,i2,i3,bsz = find3(mesh.S)
   h1,h2,h3     = mesh.h
   x1,x2,x3     = mesh.x0

   # Roman's code starts the OcTree at the top corner. Change from bottom
   # corner.
   i3 = m3 + 2 .- i3 - bsz
   x3 = x3 + m3 * h3

   S = sub2ind( (m1,m2,m3), i1,i2,i3 )
   p = sortpermFast(S)[1]
   n = length(bsz)


   # Write OcTree mesh
   f = open(fname, "w")
   println(f, m1, " ", m2, " ", m3, " ! # of cells in underlying mesh")
   println(f, x1, " ", x2, " ", x3, " ! top corner")
   println(f, h1, " ", h2, " ", h3, " ! cell size")
   println(f, n, " ! size of octree mesh")
   for i = 1:n
      idx = p[i]
      @printf(f,"%i %i %i %i\n", i1[idx], i2[idx], i3[idx], bsz[idx])
   end

   close(f)
   return
end  # function exportUBCOcTreeMesh


"""
    exportUBCOcTreeModel(fname, mesh, model)
    
    Writes an OcTree model to disk in UBC format.

    Input:
    
        name::AbstractString - File to write
        mesh::OcTreeMesh     - The mesh corresponding to the model
        model::Union{Array{Float64,1}, Array{Int64,1}} - The model 
            
"""
function exportUBCOcTreeModel(fname::AbstractString, mesh::OcTreeMesh, 
                              model::Union{Array{Float64,1}, Array{Int64,1}})
# Export OcTree cell property for use with Roman's Fortran codes.

   m1,m2,m3     = mesh.n
   i1,i2,i3,bsz = find3(mesh.S)

   # Roman's code starts the OcTree at the top corner. Change from bottom
   # corner.
   i3 = m3 + 2 .- i3 - bsz

   n = nnz(mesh.S)
   S = sub2ind( (m1,m2,m3), i1,i2,i3 )
   p = sortpermFast(S)[1]

   modelPerm = model[p]

   # Write model vector
   f = open(fname, "w")
   for i = 1:n
   	println(f, modelPerm[i])
   end
   close(f)

   return
end  # function exportUBCOcTreeModel

#-----------------------------------------------

function outputOctreeMesh(name::AbstractString, mesh::OcTreeMesh)
i1,i2,i3,bsz = find3(mesh.S)
n = nnz(mesh.S)

f = open(name, "w")
for i=1:n
	println(f, i1[i], " ", i2[i], " ", i3[i], " ", bsz[i])
end
close(f)

end  # function outputOctreeMesh



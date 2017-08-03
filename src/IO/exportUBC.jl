
export exportUBCOcTreeMesh, exportUBCOcTreeModel


"""
    exportUBCOcTreeMesh(fname, mesh)
    
    Writes an OcTree mesh to disk in UBC format.

    Input:
    
        name::AbstractString - File to write
        mesh::OcTreeMesh     - The mesh to export
            
"""
function exportUBCOcTreeMesh(fname::AbstractString, mesh::OcTreeMesh)

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
end

"""
    exportUBCOcTreeModel(fname, mesh, model)
    
    Writes an OcTree model to disk in UBC format.

    Input:
    
        name::AbstractString - File to write
        mesh::OcTreeMesh     - The mesh corresponding to the model
        model::Array{T,N}    - The model (N = 1,2)
            
"""
function exportUBCOcTreeModel{T}(name::AbstractString, mesh::OcTreeMesh, model::Array{T,1})

    m1,m2,m3     = mesh.n
    i1,i2,i3,bsz = find3(mesh.S)

    # Roman's code starts the OcTree at the top corner. Change from bottom
    # corner.
    i3 = m3 + 2 .- i3 - bsz
    n = nnz(mesh.S)
    S = sub2ind( (m1,m2,m3), i1,i2,i3 )
    p = sortpermFast(S)[1]

    # Write model vector
    f = open(name, "w")
    for i = 1:n
        idx = p[i]
        println(f, model[idx])
    end
    close(f)
    return
end

function exportUBCOcTreeModel{T}(name::AbstractString, mesh::OcTreeMesh, model::Array{T,2})

    m1,m2,m3     = mesh.n
    i1,i2,i3,bsz = find3(mesh.S)

    # Roman's code starts the OcTree at the top corner. Change from bottom
    # corner.
    i3 = m3 + 2 .- i3 - bsz
    n = nnz(mesh.S)
    S = sub2ind( (m1,m2,m3), i1,i2,i3 )
    p = sortpermFast(S)[1]
    
    ncol = size(model,2)

    # Write model vector
    f = open(name, "w")
    for i = 1:n
        idx = p[i]
        line = join(String[ string(model[idx,j]) for j = 1:ncol ], " ")
        println(f, line)
    end
    close(f)
    return
end

function outputOctreeMesh(name::AbstractString, mesh::OcTreeMesh)
    i1, i2, i3, bsz = find3(mesh.S)
    n = nnz(mesh.S)
    f = open(name, "w")
    for i = 1:n
        @printf(f,"%i %i %i %i\n", i1[i], i2[i], i3[i], bsz[i])
    end
    close(f)
end

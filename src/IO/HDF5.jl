export exportHDF5OcTreeMesh, importHDF5OcTreeMesh

"""
    exportHDF5OcTreeMesh(filename::String, mesh::OcTreeMesh; <keyword arguments>)
    exportHDF5OcTreeMesh(filename::String, meshes::Vector{T} where T<:OcTreeMesh; <keyword arguments>)
    exportHDF5OcTreeMesh(filename::String, meshes::Dict; <keyword arguments>)
    exportHDF5OcTreeMesh(filename::String, S::SparseArray3D, h::Vector{Float64}, x0::Vector{Float64}; <keyword arguments>)

Write OcTree mesh(es) to HDF5 file.

`exportHDF5OcTreeMesh(filename::String, meshes::Vector)` labels the meshes with indices
`1:length(meshes)` such that `importHDF5OcTreeMesh(filename, k)` returns `meshes[k]`.

Keyword arguments
- `compressAlg::String="blosc"`: compression algorithm; see HDF5.jl for options
- `compressLevel::Integer=4`: compression level; see HDF5.jl for options
"""

function exportHDF5OcTreeMesh(filename::String,
                              mesh::OcTreeMesh;
                              compressAlg::String="blosc",
                              compressLevel::Integer=4)

    # Create output file and mesh group
    fid = h5open(filename, "w")
    getSingleMeshHDF5group!(fid,mesh,"1",compressAlg,compressLevel)
    close(fid)
end

function exportHDF5OcTreeMesh(filename::String,
                              meshes::Vector{T} where T<:OcTreeMesh;
                              compressAlg::String="blosc",
                              compressLevel::Integer=4)
    fid = h5open(filename, "w")
    for id = 1:length(meshes)
        getSingleMeshHDF5group!(fid,meshes[id],id,compressAlg,compressLevel)
    end
    close(fid)
end

function exportHDF5OcTreeMesh(filename::String,
                              meshes::Dict;
                              compressAlg::String="blosc",
                              compressLevel::Integer=4)
    fid = h5open(filename, "w")
    for mesh in meshes
        getSingleMeshHDF5group!(fid,mesh[2],mesh[1],compressAlg,compressLevel)
    end
    close(fid)
end

function exportHDF5OcTreeMesh(filename::String,
                              S::SparseArray3D,
                              h::Vector{Float64},
                              x0::Vector{Float64};
                              compressAlg::String="blosc",
                              compressLevel::Integer=4)

    # Create output file and mesh group
    fid = h5open(filename, "w")
    getSingleMeshHDF5group!(fid,S,h,x0,"1",compressAlg,compressLevel)
    close(fid)
end

function getSingleMeshHDF5group!(fid,juliaMesh::OcTreeMesh,id,compressAlg,compressLevel)
    hdf5Mesh = g_create(fid, string(id))
    hdf5Mesh["id"] = string(id)
    attrs(hdf5Mesh)["isMesh"] = "true"

    # Add underlying mesh info
    hdf5Mesh["n"]  = collect(juliaMesh.n)
    hdf5Mesh["h"]  = juliaMesh.h
    hdf5Mesh["x0"] = juliaMesh.x0
    hdf5Mesh["nc"] = juliaMesh.nc

    # Add sparse3 as a subgroup of mesh
    s3        = g_create(hdf5Mesh,"sparse3")
    s3["i",compressAlg,compressLevel]   = juliaMesh.S.SV.nzind
    s3["bsz",compressAlg,compressLevel] = juliaMesh.S.SV.nzval
    return hdf5Mesh
end

function getSingleMeshHDF5group!(fid,S::SparseArray3D,h::Vector{Float64},x0::Vector{Float64},id,compressAlg,compressLevel)
    hdf5Mesh = g_create(fid, string(id))
    hdf5Mesh["id"] = string(id)
    attrs(hdf5Mesh)["isMesh"] = "true"

    # Add underlying mesh info
    hdf5Mesh["n"]  = [size(S)...]
    hdf5Mesh["h"]  = h
    hdf5Mesh["x0"] = x0
    hdf5Mesh["nc"] = nnz(S)

    # Add sparse3 as a subgroup of mesh
    s3        = g_create(hdf5Mesh,"sparse3")
    s3["i",compressAlg,compressLevel]   = S.SV.nzind
    s3["bsz",compressAlg,compressLevel] = S.SV.nzval
    return hdf5Mesh
end

"""
    importHDF5OcTreeMesh(filename::String)
    importHDF5OcTreeMesh(filename::String, id::N where N<:Integer)
    importHDF5OcTreeMesh(filename::String, id::String)
    importHDF5OcTreeMesh(filename::String, idList::Vector{N} where N<:Integer)
    importHDF5OcTreeMesh(filename::String, idList::Vector{String})

Import OcTree mesh(es) from HDF5 file.

`importHDF5OcTreeMesh(filename)` returns the OcTree mesh if the file contains a single
mesh. Otherwise, a dictionary with all OcTree meshes contained in the file is returned.

`importHDF5OcTreeMesh(filename, id)` returns the OcTree mesh matching the label `id`.

`importHDF5OcTreeMesh(filename, idList)` returns a dictionary with keys `idList` and
the corresponding OcTree meshes as values.
"""
function importHDF5OcTreeMesh(filename::String)
    fid = h5open(filename,"r")
    nMeshes = 0
    meshes = Dict()
    for id in names(fid)
        gmesh = fid[id]
        if exists(attrs(gmesh),"isMesh")
            nMeshes += 1
            n   = read(gmesh["n"])
            h   = read(gmesh["h"])
            x0  = read(gmesh["x0"])
            nc  = read(gmesh["nc"])
            svi = read(gmesh["sparse3/i"])
            bsz = read(gmesh["sparse3/bsz"])
            SV = sparsevec(svi,bsz,prod(n))
            S = SparseArray3D(SV,(n[1],n[2],n[3]))
            meshes[id] = getOcTreeMeshFV(S,h;x0=x0)
        else
            warn("Found non-mesh object in input file")
        end
    end
    close(fid)
    if nMeshes == 0
        error("Input file $filename did not contain any OcTree meshes")
    end
    if nMeshes == 1
        meshes = first(meshes).second
    end
    return meshes
end

importHDF5OcTreeMesh(filename::String, id::N where N<:Integer) = importHDF5OcTreeMesh(filename, string(id))

function importHDF5OcTreeMesh(filename::String, id::String)
    fid = h5open(filename,"r")
    if !any(names(fid) .== id)
      close(fid)
      error("File $filename doesn't contain mesh with id = $id")
    end
    gmesh = fid[id]
    if exists(attrs(gmesh),"isMesh")
        n   = read(gmesh["n"])
        h   = read(gmesh["h"])
        x0  = read(gmesh["x0"])
        nc  = read(gmesh["nc"])
        svi = read(gmesh["sparse3/i"])
        bsz = read(gmesh["sparse3/bsz"])
        SV = sparsevec(svi,bsz,prod(n))
        S = SparseArray3D(SV,(n[1],n[2],n[3]))
        mesh = getOcTreeMeshFV(S,h;x0=x0)
    else
        warn("Found non-mesh object in input file")
    end
    close(fid)
    return mesh
end

importHDF5OcTreeMesh(filename::String, idList::Vector{N} where N<:Integer) = importHDF5OcTreeMesh(filename, string.(idList))

function importHDF5OcTreeMesh(filename::String, idList::Vector{String})
    fid = h5open(filename,"r")
    mis = (indexin(idList, names(fid)) .== 0)
    if any(mis)
      idMis = idList[mis]
      str = idMis[1]
      for id in idMis[2:end]
        str = str * ", " * id
      end
      close(fid)
      error("File $filename doesn't contain meshes with ids = $str")
    end
    meshes = Dict()
    for id in idList
        gmesh = fid[id]
        if exists(attrs(gmesh),"isMesh")
            n   = read(gmesh["n"])
            h   = read(gmesh["h"])
            x0  = read(gmesh["x0"])
            nc  = read(gmesh["nc"])
            svi = read(gmesh["sparse3/i"])
            bsz = read(gmesh["sparse3/bsz"])
            SV = sparsevec(svi,bsz,prod(n))
            S = SparseArray3D(SV,(n[1],n[2],n[3]))
            meshes[id] = getOcTreeMeshFV(S,h;x0=x0)
        else
            warn("Found non-mesh object in input file")
        end
    end
    close(fid)
    return meshes
end

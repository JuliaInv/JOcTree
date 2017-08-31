
export getEdgeSizeNumbering, getEdgeSize, getEdgeNumbering,
       getFaceSizeNumbering, getFaceSize, getFaceNumbering,
       getNodalNumbering, getCellNumbering, getVolume, getLength,
       getVolumeVector



function getEdgeSizeNumbering(M::OcTreeMesh)
   if nnz(M.EX) != M.ne[1] || nnz(M.EY) != M.ne[2] || nnz(M.EZ) != M.ne[3] ||
      isempty(M.NEX) || isempty(M.NEY) || isempty(M.NEZ)
      M.EX,M.EY,M.EZ, M.NEX,M.NEY,M.NEZ = getEdgeSizeNumbering(M.S)
   end
   return M.EX,M.EY,M.EZ, M.NEX,M.NEY,M.NEZ
end

function getEdgeSize(M::OcTreeMesh)
   if nnz(M.EX) != M.ne[1] || nnz(M.EY) != M.ne[2] || nnz(M.EZ) != M.ne[3] ||
      isempty(M.NEX) || isempty(M.NEY) || isempty(M.NEZ)
      M.EX,M.EY,M.EZ, M.NEX,M.NEY,M.NEZ = getEdgeSizeNumbering(M.S)
   end
   return M.EX,M.EY,M.EZ
end

function getEdgeSize(S::SparseArray3D)
   EX,EY,EZ, NEX,NEY,NEZ = getEdgeSizeNumbering(S)
   return EX,EY,EZ
end

function getEdgeNumbering(M::OcTreeMesh)
   if nnz(M.EX) != M.ne[1] || nnz(M.EY) != M.ne[2] || nnz(M.EZ) != M.ne[3] ||
      isempty(M.NEX) || isempty(M.NEY) || isempty(M.NEZ)
      M.EX,M.EY,M.EZ, M.NEX,M.NEY,M.NEZ = getEdgeSizeNumbering(M.S)
   end
   return M.NEX,M.NEY,M.NEZ
end

function getEdgeNumbering(S::SparseArray3D)
   EX,EY,EZ, NEX,NEY,NEZ = getEdgeSizeNumbering(S)
   return NEX,NEY,NEZ
end


function getEdgeSizeNumbering(S::SparseArray3D)

 m1,m2,m3 = S.sz
 i,j,k,bsz = find3(S)

 ns = nnz(S)
 ns4 = ns*4
 ii = Array{Int64}(ns4)
 jj = Array{Int64}(ns4)
 kk = Array{Int64}(ns4)
 vv = Array{Int64}(ns4)

 vv[1:4:ns4] = bsz
 vv[2:4:ns4] = bsz
 vv[3:4:ns4] = bsz
 vv[4:4:ns4] = bsz

 # X edges
 sizeEX      = (m1,m2+1,m3+1)
 ii[1:4:ns4] = i
 ii[2:4:ns4] = i
 ii[3:4:ns4] = i
 ii[4:4:ns4] = i

 jj[1:4:ns4] = j
 jj[2:4:ns4] = j + bsz
 jj[3:4:ns4] = j
 jj[4:4:ns4] = j + bsz

 kk[1:4:ns4] = k
 kk[2:4:ns4] = k
 kk[3:4:ns4] = k + bsz
 kk[4:4:ns4] = k + bsz

 EX  = sparse3(ii,jj,kk, vv, collect(sizeEX))
 EXN = deepcopy(EX)
 copy!(EXN.SV.nzval, 1:nnz(EX) )


 # Y edges
 sizeEY      = (m1+1,m2,m3+1)
 ii[1:4:ns4] = i
 ii[2:4:ns4] = i + bsz
 ii[3:4:ns4] = i
 ii[4:4:ns4] = i + bsz

 jj[1:4:ns4] = j
 jj[2:4:ns4] = j
 jj[3:4:ns4] = j
 jj[4:4:ns4] = j

 kk[1:4:ns4] = k
 kk[2:4:ns4] = k
 kk[3:4:ns4] = k + bsz
 kk[4:4:ns4] = k + bsz

 EY  = sparse3(ii,jj,kk, vv, collect(sizeEY))
 EYN = deepcopy(EY)
 copy!(EYN.SV.nzval, 1:nnz(EY) )


 # Z edges
 sizeEZ      = (m1+1,m2+1,m3)
 ii[1:4:ns4] = i
 ii[2:4:ns4] = i + bsz
 ii[3:4:ns4] = i
 ii[4:4:ns4] = i + bsz

 jj[1:4:ns4] = j
 jj[2:4:ns4] = j
 jj[3:4:ns4] = j + bsz
 jj[4:4:ns4] = j + bsz

 kk[1:4:ns4] = k
 kk[2:4:ns4] = k
 kk[3:4:ns4] = k
 kk[4:4:ns4] = k

 EZ  = sparse3(ii,jj,kk, vv, collect(sizeEZ))
 EZN = deepcopy(EZ)
 copy!(EZN.SV.nzval, 1:nnz(EZ) )


return EX,  EY,  EZ,   # edge sizes
       EXN, EYN, EZN   # edge numbering
end  # function getEdgeSizeNumbering

"""
Mesh.L = getLength(Mesh::OcTreeMesh) computes edge lengths l, returns spdiagm(l)
"""
function getLength(Mesh::OcTreeMesh)
    if isempty(Mesh.L)
        l1, l2, l3 = getEdgeSize(Mesh)
        l1 = nonzeros(l1)*Mesh.h[1]
        l2 = nonzeros(l2)*Mesh.h[2]
        l3 = nonzeros(l3)*Mesh.h[3]
        Mesh.L = spdiagm([l1;l2;l3])
    end
    return Mesh.L
end

#-------------------------------------------------------------------------

function getFaceSizeNumbering(M::OcTreeMesh)
   if nnz(M.FX) != M.nf[1] || nnz(M.FY) != M.nf[2] || nnz(M.FZ) != M.nf[3] ||
      isempty(M.NFX) || isempty(M.NFY) || isempty(M.NFZ)
      M.FX,M.FY,M.FZ, M.NFX,M.NFY,M.NFZ  = getFaceSizeNumbering(M.S)
   end
   return M.FX,M.FY,M.FZ, M.NFX,M.NFY,M.NFZ
end

function getFaceSize(M::OcTreeMesh)
   if nnz(M.FX) != M.nf[1] || nnz(M.FY) != M.nf[2] || nnz(M.FZ) != M.nf[3] ||
      isempty(M.NFX) || isempty(M.NFY) || isempty(M.NFZ)
      M.FX,M.FY,M.FZ, M.NFX,M.NFY,M.NFZ = getFaceSizeNumbering(M.S)
   end
   return M.FX,M.FY,M.FZ
end

function getFaceSize(S::SparseArray3D)
   FX,FY,FZ, NFX,NFY,NFZ = getFaceSizeNumbering(S)
   return FX,FY,FZ
end

function getFaceNumbering(M::OcTreeMesh)
   if nnz(M.FX) != M.nf[1] || nnz(M.FY) != M.nf[2] || nnz(M.FZ) != M.nf[3] ||
      isempty(M.NFX) || isempty(M.NFY) || isempty(M.NFZ)
      M.FX,M.FY,M.FZ, M.NFX,M.NFY,M.NFZ  = getFaceSizeNumbering(M.S)
   end
   return M.NFX,M.NFY,M.NFZ
end

function getFaceNumbering(S::SparseArray3D)
   FX,FY,FZ, NFX,NFY,NFZ = getFaceSizeNumbering(S)
   return NFX,NFY,NFZ
end


function getFaceSizeNumbering(S::SparseArray3D)

 m1,m2,m3 = S.sz
 i,j,k,bsz = find3(S)

 ns = nnz(S)
 ns2 = ns*2
 ii = Array{Int64}(ns2)
 jj = Array{Int64}(ns2)
 kk = Array{Int64}(ns2)
 vv = Array{Int64}(ns2)

 vv[1:2:ns2] = bsz
 vv[2:2:ns2] = bsz


 # X faces
 sizeFX = (m1+1, m2, m3)

 ii[1:2:ns2] = i
 ii[2:2:ns2] = i + bsz
 jj[1:2:ns2] = j
 jj[2:2:ns2] = j
 kk[1:2:ns2] = k
 kk[2:2:ns2] = k

 FX  = sparse3(ii,jj,kk, vv, collect(sizeFX))
 FXN = deepcopy(FX)
 copy!(FXN.SV.nzval, 1:nnz(FX) )


 # Y faces
 sizeFY = (m1, m2+1, m3)

 ii[1:2:ns2] = i
 ii[2:2:ns2] = i
 jj[1:2:ns2] = j
 jj[2:2:ns2] = j + bsz
 kk[1:2:ns2] = k
 kk[2:2:ns2] = k

 FY  = sparse3(ii,jj,kk, vv, collect(sizeFY))
 FYN = deepcopy(FY)
 copy!(FYN.SV.nzval, 1:nnz(FY) )


 # Z faces
 sizeFZ = (m1, m2, m3+1)
 ii[1:2:ns2] = i
 ii[2:2:ns2] = i
 jj[1:2:ns2] = j
 jj[2:2:ns2] = j
 kk[1:2:ns2] = k
 kk[2:2:ns2] = k + bsz

 FZ  = sparse3(ii,jj,kk, vv, collect(sizeFZ))
 FZN = deepcopy(FZ)
 copy!(FZN.SV.nzval, 1:nnz(FZ) )


return FX,  FY,  FZ,   # face sizes
       FXN, FYN, FZN   # face numbering
end  # function getFaceSizeNumbering

#-------------------------------------------------------------------------

function getNodalNumbering(M::OcTreeMesh)
   if (nnz(M.NN) != M.nn) || isempty(M.NN)
      M.NN = getNodalNumbering(M.S)
   end
   return M.NN
end

function getNodalNumbering(S::SparseArray3D)
    # Numbering of the nodes of an OcTree structure

    m1,m2,m3 = S.sz
    i,j,k,bsz = find3(S)

    ns = nnz(S)
    ns8 = ns*8
    ii = Array{Int64}(ns8)
    jj = Array{Int64}(ns8)
    kk = Array{Int64}(ns8)

    ii[1:8:ns8] = i
    ii[2:8:ns8] = i + bsz
    ii[3:8:ns8] = i
    ii[4:8:ns8] = i + bsz
    ii[5:8:ns8] = i
    ii[6:8:ns8] = i + bsz
    ii[7:8:ns8] = i
    ii[8:8:ns8] = i + bsz

    jj[1:8:ns8] = j
    jj[2:8:ns8] = j
    jj[3:8:ns8] = j + bsz
    jj[4:8:ns8] = j + bsz
    jj[5:8:ns8] = j
    jj[6:8:ns8] = j
    jj[7:8:ns8] = j + bsz
    jj[8:8:ns8] = j + bsz

    kk[1:8:ns8] = k
    kk[2:8:ns8] = k
    kk[3:8:ns8] = k
    kk[4:8:ns8] = k
    kk[5:8:ns8] = k + bsz
    kk[6:8:ns8] = k + bsz
    kk[7:8:ns8] = k + bsz
    kk[8:8:ns8] = k + bsz

    N = sparse3(ii,jj,kk, kk, [m1+1,m2+1,m3+1])
    copy!(N.SV.nzval, 1:nnz(N) )

    return N
end  # function getNodalNumbering

#-------------------------------------------------------------------------

function getCellNumbering(M::OcTreeMesh)
  if (nnz(M.NC) != M.nc) || isempty(M.NC)
     M.NC = getCellNumbering(M.S)
  end
  return M.NC
end

function getCellNumbering(S::SparseArray3D)
    NC = deepcopy(S)
    copy!(NC.SV.nzval, 1:nnz(S) )
    return NC
end  # function getCellNumbering

"""
Mesh.V = getVolume(M::OcTreeMesh) returns diagonal matrix of cell volumes
"""
function getVolume(M::OcTreeMesh)
    if isempty(M.V)
        i,j,k,bsz = find3(M.S)
        h         = M.h
        M.V       = spdiagm(bsz.^3*prod(h))
    end
    return M.V
end

function getVolumeVector(M::OcTreeMesh)
    if isempty(M.V)
        i,j,k,bsz = find3(M.S)
        h         = M.h
        M.V       = spdiagm(bsz.^3*prod(h))
    end
    return M.V.nzval
end

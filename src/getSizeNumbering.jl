
export getEdgeSizeNumbering, getEdgeSize, getEdgeNumbering,
       getFaceSizeNumbering, getFaceSize, getFaceNumbering



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


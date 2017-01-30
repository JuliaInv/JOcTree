export getCellNumbering

function getCellNumbering(S::SparseArray3D)

sz = collect(size(S.SV))
colptr = copy(S.SV.colptr)
rowval = copy(S.SV.rowval)

nz = length(rowval)
nzval = collect(1:nz)

SV = SparseMatrixCSC(sz[1],sz[2], colptr, rowval, nzval)

CN = SparseArray3D(SV, S.sz)
return CN
end  # function getCellNumbering

export getBoundaryEdges

function getBoundaryEdges(S::SparseArray3D)
# [BEX,BEY,BEZ,BE] = getBoundaryEdges(S)
#
# Purpose: Find the edges which are on the boundary
#
# ----------------------------------------------------------------------------
 
EXN, EYN, EZN = getEdgeNumbering(S)
 
#BEX = uniqueidx([ nonzeros(EXN[:,:,[1 end]]); nonzeros(EXN[:,[1 end],:]) ])[1]
s1,s2,s3 = size(EXN)
i,j,k = find3(EXN)
BEX = [ find(j.==1) ; find(j.==s2) ; find(k.==1) ; find(k.==s3) ]


#BEY = uniqueidx([ nonzeros(EYN[[1 end],:,:]); nonzeros(EYN[:,:,[1 end]]) ])[1]
s1,s2,s3 = size(EYN)
i,j,k = find3(EYN)
BEY = [ find(i.==1) ; find(i.==s1) ; find(k.==1) ; find(k.==s3) ]


#BEZ = uniqueidx([ nonzeros(EZN[[1 end],:,:]); nonzeros(EZN[:,[1 end],:]) ])[1]
s1,s2,s3 = size(EZN)
i,j,k = find3(EZN)
BEZ = [ find(i.==1) ; find(i.==s1) ; find(j.==1) ; find(j.==s2) ]
 
 
BE = [ BEX ; BEY + nnz(EXN) ; BEZ + nnz(EXN) + nnz(EYN) ]

return BEX,BEY,BEZ, BE
end  # function getBoundaryEdges

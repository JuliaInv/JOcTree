
# NOT USED

export getCellNumbering  #, getNodalNumbering # getEdgeNumbering, getFaceNumbering, 

function getCellNumbering(M::OcTreeMesh)
    return getCellNumbering(M.S)
end

function getCellNumbering(S::SparseArray3D)
    sz = collect(size(S.SV))
    colptr = copy(S.SV.colptr)
    rowval = copy(S.SV.rowval)

    nz = length(rowval)
    nzval = collect(1:nz)

    SV = SparseMatrixCSC(sz[1],sz[2], colptr, rowval, nzval)

    CN = SparseArray3D(SV, S.sz)
    return CN
end

# ------------------------------------------------------------

#function getEdgeNumbering(M::OcTreeMesh)
#    if isempty(M.NEX) || isempty(M.NEY) || isempty(M.NEZ)
#        M.NEX,M.NEY,M.NEZ = getEdgeNumbering(M.S)
#    end
#    return M.NEX,M.NEY,M.NEZ
#end
#
#function getEdgeNumbering(S)
#    # [EX,EY,EZ] = getEdgeNumbering(S)
#    # Numbering of the edges of the structure
#
#    m1,m2,m3 = S.sz
#    i,j,k,bsz = find3(S)
#
#
#    sizeEX      = (m1,m2+1,m3+1)
#    ii          = [  i      ;   i       ;   i       ;   i       ;]
#    jj          = [  j      ;   j+bsz   ;   j       ;   j+bsz   ;]
#    kk          = [  k      ;   k       ;   k+bsz   ;   k+bsz   ;]
#    ii,jj,kk  = ind2sub(sizeEX,sort(unique(sub2ind(sizeEX,ii,jj,kk))))  # make em unique
#    EX          = sparse3(ii,jj,kk,1:length(ii), collect(sizeEX))
#
#    sizeEY      = (m1+1,m2,m3+1)
#    ii          = [  i      ;   i+bsz   ;   i       ;   i+bsz   ;]
#    jj          = [  j      ;   j       ;   j       ;   j       ;]
#    kk          = [  k      ;   k       ;   k+bsz   ;   k+bsz   ;]
#    ii,jj,kk  = ind2sub(sizeEY,sort(unique(sub2ind(sizeEY,ii,jj,kk))))  # make em unique
#    EY          = sparse3(ii,jj,kk,1:length(ii), collect(sizeEY))
#
#    sizeEZ      = (m1+1,m2+1,m3)
#    ii          = [  i       ;    i+bsz   ;   i       ;   i+bsz   ;]
#    jj          = [  j       ;    j       ;   j+bsz   ;   j+bsz   ;]
#    kk          = [  k       ;    k       ;   k       ;   k       ;]
#    ii,jj,kk  = ind2sub(sizeEZ,sort(unique(sub2ind(sizeEZ,ii,jj,kk))))  # make em unique
#    EZ          = sparse3(ii,jj,kk,1:length(ii), collect(sizeEZ))
#
#    return EX, EY, EZ
#end

# ------------------------------------------------------------

#function getFaceNumbering(M::OcTreeMesh)
#    if isempty(M.NFX) || isempty(M.NFY) || isempty(M.NFZ)
#        M.NFX,M.NFY,M.NFZ = getFaceNumbering(M.S)
#    end
#    return M.NFX,M.NFY,M.NFZ
#end
#
#function getFaceNumbering(S)
#
#    m1,m2,m3 = S.sz
#    i,j,k,bsz = find3(S)
#
#    ##     mark upper    mark lower
#    ##      |                |
#    ##      v                v
#    ii = vcat(  i          , i+bsz )
#    jj = vcat(  j          , j     )
#    kk = vcat(  k          , k     )
#
#    sizeFX = (m1+1, m2, m3)
#    I = sort(unique(sub2ind(sizeFX, ii,jj,kk)))   ## create unique sorted linear indices
#    ii,jj,kk = ind2sub(sizeFX, I)         ## linear indices to nd indices
#    FX = sparse3(ii,jj,kk,1:length(ii), collect(sizeFX))
#
#    ##     mark left     mark right
#    ##      |                |
#    ##      v                v
#    ii = vcat(  i          , i    )
#    jj = vcat(  j          , j+bsz)
#    kk = vcat(  k          , k    )
#    sizeFY = (m1, m2+1, m3)
#    I = sort(unique(sub2ind(sizeFY, ii,jj,kk)))   ## create unique sorted linear indices
#    ii,jj,kk = ind2sub(sizeFY, I)         ## linear indices to nd indices
#    FY = sparse3(ii,jj,kk,1:length(ii), collect(sizeFY))
#
#    ##     mark front    mark back
#    ##      |                |
#    ##      v                v
#    ii = vcat(  i          , i     )
#    jj = vcat(  j          , j     )
#    kk = vcat(  k          , k+bsz )
#
#    sizeFZ = (m1, m2, m3+1)
#    I = sort(unique(sub2ind(sizeFZ, ii,jj,kk)))   ## create unique sorted linear indices
#    ii,jj,kk = ind2sub(sizeFZ, I)         ## linear indices to nd indices
#    FZ = sparse3(ii,jj,kk,1:length(ii), collect(sizeFZ))
#
#    return FX, FY, FZ
#end

# ------------------------------------------------------------

#function getNodalNumbering(M::OcTreeMesh)
#    return getNodalNumbering(M.S)
#end
#
#function getNodalNumbering(S)
#    # N = getEdgeNumbering(S)
#    # Numbering of the nodes of an OcTree structure
#
#    m1,m2,m3 = S.sz
#    i,j,k,bsz = find3(S)
#
#    nind = [i        j       k;
#            i        j       k+bsz;
#            i        j+bsz   k;
#            i        j+bsz   k+bsz;
#            i+bsz    j       k;
#            i+bsz    j       k+bsz;
#            i+bsz    j+bsz   k;
#            i+bsz    j+bsz   k+bsz ]
#
#    Ni = sparse3(nind[:,1],nind[:,2],nind[:,3],nind[:,3],[m1+1,m2+1,m3+1])
#
#    i,j,k = find3(Ni)
#    N = sparse3(i,j,k,1:length(i), [m1+1,m2+1,m3+1]);
#
#    return N
#end

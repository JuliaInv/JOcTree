export regularizeOcTree


function regularizeOcTree(S::SparseArray3D)
# S = regularizeOcTree(S)

isRegularTree = false

while !isRegularTree
    isRegularTree = true

    #i1,j1,k1,bsz1,i2,j2,k2,bsz2 = findNonRegularBlocks_(S)
    Inr, i,j,k, bsz = findNonRegularBlocks(S)

    if !isempty(Inr)
       # compute entries that stays unmodified (Ic = set complement of I)
       isRegularTree = false

        #bsz1 = div(bsz1, 2)
        

         #m1,m2,m3 = S.sz
         
         #ii = vcat(i1, i1+bsz1, i1,      i1+bsz1, i1,      i1+bsz1, i1,      i1+bsz1, i2)
         #jj = vcat(j1, j1,      j1+bsz1, j1+bsz1, j1,      j1,      j1+bsz1, j1+bsz1, j2)
         #kk = vcat(k1, k1,      k1,      k1,      k1+bsz1, k1+bsz1, k1+bsz1, k1+bsz1, k2)
         #vv = vcat(bsz1, bsz1,  bsz1,    bsz1,    bsz1,    bsz1,    bsz1,    bsz1,   bsz2)
         

         #S = sparse3(ii,jj,kk,vv,[m1,m2,m3])
         
         S = splitCells(i,j,k,bsz, S.sz, Inr)
    end

end

return S

end # function regularizeOcTree

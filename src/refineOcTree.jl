export refineOcTree

function refineOcTree(S::SparseArray3D, tau::Vector, tol)

i,j,k,bsz = find3(S)


#I = find([bsz .>1] .== true & [abs(tau) .>tol] .== true)
#I = find([bsz .>1;] .== true & [abs(tau) .>tol;] .== true)
I = find( (bsz .> 1)  & (abs(tau) .> tol) )

if !isempty(I)
    # compute entries stying unmodified (Ic = set complement of I)
    Ic = setdiff(vec(1:length(bsz)) , vec(I))

    bsz[I] = div(bsz[I], 2)


    m1,m2,m3 = S.sz

    # modified/new entries  followed by  unmodified old entries
    ii = vcat(i[I], i[I]+bsz[I], i[I],        i[I]+bsz[I], i[I],        i[I]+bsz[I], i[I],        i[I]+bsz[I], i[Ic])
    jj = vcat(j[I], j[I],        j[I]+bsz[I], j[I]+bsz[I], j[I],        j[I],        j[I]+bsz[I], j[I]+bsz[I], j[Ic])
    kk = vcat(k[I], k[I],        k[I],        k[I],        k[I]+bsz[I], k[I]+bsz[I], k[I]+bsz[I], k[I]+bsz[I], k[Ic])
    vv = vcat(bsz[I], bsz[I], bsz[I], bsz[I], bsz[I], bsz[I], bsz[I], bsz[I], bsz[Ic])

    SF = sparse3(ii,jj,kk,vv,[m1,m2,m3])
else
    SF = S
end

return SF

end  # function refineOcTree

export getNodalConstraints

function getNodalConstraints(M::OcTreeMesh)
    if isempty(M.Nn)
        T = eltype(M.h)
        if all(M.S.SV.nzval.==M.S.SV.nzval[1]) # uniform mesh
            N = typeof(M.nc)
            Nn = prod(M.n+1)
            M.Nn = speye(Nn)
            M.Qn = speye(Nn)
        else
            M.Nn,M.Qn,Cn,M.activeNodes = getNodalConstraints(T,M.S)
        end
    end
    return M.Nn,M.Qn,M.activeNodes
end

function getNodalConstraints(::Type{Tf},S::SparseArray3D) where Tf
    i,j,k,bsz = find3(S)


    upper,lower,left,right,front,back = getNumberOfNeighbors(S)
    nn = [upper; lower; left; right; front; back]
    if ~all( (nn .== 0) .| (nn .== 1) .| (nn .== 4) )
    error("Implemented only for regularized OcTree meshes")
    end
    NN = getNodalNumbering(S)
    Tn  = eltype(NN.SV.nzval)
    Tn2 = eltype(NN.SV.nzind)
    bsz = convert(Vector{Tn2},bsz)
    n1 = Tn[]; n2 = Tn[]; n3 = Tn[];
    n4 = Tn[]; n5 = Tn[]; n6 = Tn[];
    n7 = Tn[]; n8 = Tn[]; n9 = Tn[];

    #    ^ k z
    #   /
    #  /
    # +-------> j y
    # |
    # |
    # v x


    #         n7--------n8--------n9
    #         /         /         /
    #        /         /         /
    #       /         /         /
    #     n4--------n5--------n6
    #     /         /         /
    #    /         /         /
    #   /         /         /
    # n1--------n2--------n3

    I = find(upper.==4) #find "bigger" cells
    bsz2 = div.(bsz[I],Tn2(2))

    if !isempty(I)
    append!(n1, NN.SV[sub2ind(size(NN), i[I], j[I]        , k[I]        )])
    append!(n2, NN.SV[sub2ind(size(NN), i[I], j[I]+bsz2   , k[I]        )])
    append!(n3, NN.SV[sub2ind(size(NN), i[I], j[I]+bsz[I] , k[I]        )])
    append!(n4, NN.SV[sub2ind(size(NN), i[I], j[I]        , k[I]+bsz2   )])
    append!(n5, NN.SV[sub2ind(size(NN), i[I], j[I]+bsz2   , k[I]+bsz2   )])
    append!(n6, NN.SV[sub2ind(size(NN), i[I], j[I]+bsz[I] , k[I]+bsz2   )])
    append!(n7, NN.SV[sub2ind(size(NN), i[I], j[I]        , k[I]+bsz[I] )])
    append!(n8, NN.SV[sub2ind(size(NN), i[I], j[I]+bsz2   , k[I]+bsz[I] )])
    append!(n9, NN.SV[sub2ind(size(NN), i[I], j[I]+bsz[I] , k[I]+bsz[I] )])
    end

    I = find(lower.==4) # find  "bigger" cells
    bsz2 = div.(bsz[I],Tn2(2))

    if !isempty(I)
    append!(n1, NN.SV[sub2ind(size(NN), i[I]+bsz[I], j[I]        , k[I]       )])
    append!(n2, NN.SV[sub2ind(size(NN), i[I]+bsz[I], j[I]+bsz2   , k[I]       )])
    append!(n3, NN.SV[sub2ind(size(NN), i[I]+bsz[I], j[I]+bsz[I] , k[I]       )])
    append!(n4, NN.SV[sub2ind(size(NN), i[I]+bsz[I], j[I]        , k[I]+bsz2  )])
    append!(n5, NN.SV[sub2ind(size(NN), i[I]+bsz[I], j[I]+bsz2   , k[I]+bsz2  )])
    append!(n6, NN.SV[sub2ind(size(NN), i[I]+bsz[I], j[I]+bsz[I] , k[I]+bsz2  )])
    append!(n7, NN.SV[sub2ind(size(NN), i[I]+bsz[I], j[I]        , k[I]+bsz[I])])
    append!(n8, NN.SV[sub2ind(size(NN), i[I]+bsz[I], j[I]+bsz2   , k[I]+bsz[I])])
    append!(n9, NN.SV[sub2ind(size(NN), i[I]+bsz[I], j[I]+bsz[I] , k[I]+bsz[I])])
    end

    ##################################

    I = find(left .== 4) # find  "bigger" cells
    bsz2 = div.(bsz[I],Tn2(2))

    if !isempty(I)
    append!(n1, NN.SV[sub2ind(size(NN), i[I]        , j[I], k[I]        )])
    append!(n2, NN.SV[sub2ind(size(NN), i[I]+bsz2   , j[I], k[I]        )])
    append!(n3, NN.SV[sub2ind(size(NN), i[I]+bsz[I] , j[I], k[I]        )])
    append!(n4, NN.SV[sub2ind(size(NN), i[I]        , j[I], k[I]+bsz2   )])
    append!(n5, NN.SV[sub2ind(size(NN), i[I]+bsz2   , j[I], k[I]+bsz2   )])
    append!(n6, NN.SV[sub2ind(size(NN), i[I]+bsz[I] , j[I], k[I]+bsz2   )])
    append!(n7, NN.SV[sub2ind(size(NN), i[I]        , j[I], k[I]+bsz[I] )])
    append!(n8, NN.SV[sub2ind(size(NN), i[I]+bsz2   , j[I], k[I]+bsz[I] )])
    append!(n9, NN.SV[sub2ind(size(NN), i[I]+bsz[I] , j[I], k[I]+bsz[I] )])
    end

    I = find(right .== 4) # find  "bigger" cells
    bsz2 = div.(bsz[I],Tn2(2))

    if !isempty(I)
    append!(n1, NN.SV[sub2ind(size(NN), i[I]         , j[I]+bsz[I], k[I]        )])
    append!(n2, NN.SV[sub2ind(size(NN), i[I]+bsz2    , j[I]+bsz[I], k[I]        )])
    append!(n3, NN.SV[sub2ind(size(NN), i[I]+bsz[I]  , j[I]+bsz[I], k[I]        )])
    append!(n4, NN.SV[sub2ind(size(NN), i[I]         , j[I]+bsz[I], k[I]+bsz2   )])
    append!(n5, NN.SV[sub2ind(size(NN), i[I]+bsz2    , j[I]+bsz[I], k[I]+bsz2   )])
    append!(n6, NN.SV[sub2ind(size(NN), i[I]+bsz[I]  , j[I]+bsz[I], k[I]+bsz2   )])
    append!(n7, NN.SV[sub2ind(size(NN), i[I]         , j[I]+bsz[I], k[I]+bsz[I] )])
    append!(n8, NN.SV[sub2ind(size(NN), i[I]+bsz2    , j[I]+bsz[I], k[I]+bsz[I] )])
    append!(n9, NN.SV[sub2ind(size(NN), i[I]+bsz[I]  , j[I]+bsz[I], k[I]+bsz[I] )])
    end

    ##################################

    I = find(front .== 4) # find  "bigger" cells
    bsz2 = div.(bsz[I],Tn2(2))

    if !isempty(I)
    append!(n1, NN.SV[sub2ind(size(NN), i[I]         , j[I]        , k[I] )])
    append!(n2, NN.SV[sub2ind(size(NN), i[I]+bsz2    , j[I]        , k[I] )])
    append!(n3, NN.SV[sub2ind(size(NN), i[I]+bsz[I]  , j[I]        , k[I] )])
    append!(n4, NN.SV[sub2ind(size(NN), i[I]         , j[I]+bsz2   , k[I] )])
    append!(n5, NN.SV[sub2ind(size(NN), i[I]+bsz2    , j[I]+bsz2   , k[I] )])
    append!(n6, NN.SV[sub2ind(size(NN), i[I]+bsz[I]  , j[I]+bsz2   , k[I] )])
    append!(n7, NN.SV[sub2ind(size(NN), i[I]         , j[I]+bsz[I] , k[I] )])
    append!(n8, NN.SV[sub2ind(size(NN), i[I]+bsz2    , j[I]+bsz[I] , k[I] )])
    append!(n9, NN.SV[sub2ind(size(NN), i[I]+bsz[I]  , j[I]+bsz[I] , k[I] )])
    end

    I = find(back .== 4) # find  "bigger" cells
    bsz2 = div.(bsz[I],Tn2(2))

    if !isempty(I)
    append!(n1, NN.SV[sub2ind(size(NN), i[I]        , j[I]        , k[I]+bsz[I] )])
    append!(n2, NN.SV[sub2ind(size(NN), i[I]+bsz2   , j[I]        , k[I]+bsz[I] )])
    append!(n3, NN.SV[sub2ind(size(NN), i[I]+bsz[I] , j[I]        , k[I]+bsz[I] )])
    append!(n4, NN.SV[sub2ind(size(NN), i[I]        , j[I]+bsz2   , k[I]+bsz[I] )])
    append!(n5, NN.SV[sub2ind(size(NN), i[I]+bsz2   , j[I]+bsz2   , k[I]+bsz[I] )])
    append!(n6, NN.SV[sub2ind(size(NN), i[I]+bsz[I] , j[I]+bsz2   , k[I]+bsz[I] )])
    append!(n7, NN.SV[sub2ind(size(NN), i[I]        , j[I]+bsz[I] , k[I]+bsz[I] )])
    append!(n8, NN.SV[sub2ind(size(NN), i[I]+bsz2   , j[I]+bsz[I] , k[I]+bsz[I] )])
    append!(n9, NN.SV[sub2ind(size(NN), i[I]+bsz[I] , j[I]+bsz[I] , k[I]+bsz[I] )])
    end

    ###################################################################
    ##
    ## Build interpolation matrix
    ##
    ###################################################################
    # Constraint matrix
    #
    #   n1  n3  n7  n9  n2  n4  n6  n8  n5
    #  I/2 I/2          -I
    #  I/2     I/2          -I
    #      I/2     I/2          -I
    #          I/2 I/2              -I
    #  I/4 I/4 I/4 I/4                  -I
    #
    # Null space matrix
    #
    #  n1 I
    #  n3     I
    #  n7         I
    #  n9             I
    #  n2 I/2 I/2
    #  n4 I/2     I/2
    #  n6     I/2     I/2
    #  n8         I/2 I/2
    #  n5 I/4 I/4 I/4 I/4

    n = Tn(nnz(NN))
    I = speye(Tf,Tn,n)

    # eliminate nodes at edge centers (n2,n4,n6,n8)
    i1 = [n1;n1;n3;n7]
    i2 = [n3;n7;n9;n9]
    i3 = [n2;n4;n6;n8]
    i3orig = copy(i3)
    i3 = unique(i3)
    i = indexin(i3,i3orig)
    i3orig = 0.0
    i1 = i1[i]
    i2 = i2[i]
    A1 = I[i1,:]/2 + I[i2,:]/2
    B1 = I[i3,:]

    # eliminate nodes at face centers (n5)
    j  = sortperm(vec(n5))
    j5 = n5[j]
    j1 = n1[j]
    j2 = n3[j]
    j3 = n7[j]
    j4 = n9[j]
    A2 = I[j1,:]/4 + I[j2,:]/4 + I[j3,:]/4 + I[j4,:]/4;
    B2 = I[j5,:]

    A = [A1; A2]
    B = [B1; B2]

    # permute rows and columns for matrix splitting
    Tn1 = one(Tn)
    k2 = sort([i3; j5])
    p  = sortperm([i3; j5])
    k1 = setdiff([Tn1:n;], k2)
    q  = [k1; k2]
    A  = A[p,q]
    B  = B[p,q]
    m  = n - length(p);

    C = A - B;
    N = [
    speye(m)
    A[:,1:m]]
    Q = [speye(m) spzeros(m,length(p))]

    # permute back
    p = sortperm(q)
    C = C[:,p]
    N = N[p,:]
    Q = Q[:,p]

    # table lookup for new nodal numbering
    p = zeros(Tn,n)
    p[k1] = 1:length(k1)

    return N,Q,C,p

end

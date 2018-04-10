export getFaceConstraints

function getFaceConstraints(M::OcTreeMesh)
	if isempty(M.Nf)
        T = eltype(M.h)
		if all(M.S.SV.nzval.==M.S.SV.nzval[1]) # uniform mesh
			N    = typeof(M.nc)
			M.Nf = speye(T,N,sum(M.nf))
			M.Qf = speye(T,N,sum(M.nf))
			M.activeFaces = collect(N,1:sum(M.nf))
		else
			M.Nf,M.Qf,Cf,M.activeFaces = getFaceConstraints(T,M.S)
		end
	end
	return M.Nf,M.Qf,M.activeFaces
end

function getFaceConstraints(::Type{Tf},S::SparseArray3D) where Tf
    #Eliminate hanging faces
    #   C ... constraint matrix
    #   N ... null space matrix (interpolation: coarse to fine)
    #   Q ... projection matrix (restriction: fine to coarse)
    #                     which satisfy Q * N = I
    #                     (note that Q ~= N')
    #   p ............... lookup table for new face enumeration:
    #                     new face number = p(old face number)
    #                     (returns zero for deleted faces)
    i0,j0,k0,bsz = find3(S)
	Tn2 = eltype(i0)
	bsz = convert(Vector{Tn2},bsz)

    i1 = i0 + div.(bsz, Tn2(2))
    i2 = i0 + bsz
    j1 = j0 + div.(bsz, Tn2(2))
    j2 = j0 + bsz
    k1 = k0 + div.(bsz, Tn2(2))
    k2 = k0 + bsz

    upper,lower,left,right,front,back = getNumberOfNeighbors(S);
    nn = [upper; lower; left; right; front; back];
    if ~all( (nn .== 0) .| (nn .== 1) .| (nn .== 4) )
    error("Implemented only for regularized OcTree meshes")
    end

    FNX,FNY,FNZ = getFaceNumbering(S);

	Tn  = eltype(FNX.SV.nzval)
    fx1 = Tn[]; fx2 = Tn[]; fx3 = Tn[]; fx4 = Tn[];
    fy1 = Tn[]; fy2 = Tn[]; fy3 = Tn[]; fy4 = Tn[];
    fz1 = Tn[]; fz2 = Tn[]; fz3 = Tn[]; fz4 = Tn[];


    #    ^ k z
    #   /
    #  /
    # +-------> j y
    # |
    # |
    # v x
    #
    #
    # x-faces
    #         +---------+---------+
    #        /         /         /
    #       /   fx3   /   fx4   /
    #      /         /         /
    #     +---------+---------+
    #    /         /         /
    #   /   fx1   /   fx2   /
    #  /         /         /
    # +---------+---------+

    I = (upper .== 4) # find  "bigger" cells

    if any(I)
    append!(fx1, FNX.SV[sub2ind(size(FNX), i0[I], j0[I], k0[I])])
    append!(fx2, FNX.SV[sub2ind(size(FNX), i0[I], j1[I], k0[I])])
    append!(fx3, FNX.SV[sub2ind(size(FNX), i0[I], j0[I], k1[I])])
    append!(fx4, FNX.SV[sub2ind(size(FNX), i0[I], j1[I], k1[I])])
    end

    I = (lower .== 4) # find  "bigger" cells

    if any(I)
    append!(fx1, FNX.SV[sub2ind(size(FNX), i2[I], j0[I], k0[I])])
    append!(fx2, FNX.SV[sub2ind(size(FNX), i2[I], j1[I], k0[I])])
    append!(fx3, FNX.SV[sub2ind(size(FNX), i2[I], j0[I], k1[I])])
    append!(fx4, FNX.SV[sub2ind(size(FNX), i2[I], j1[I], k1[I])])
    end


    ################################

    I = (left .== 4) # find  "bigger" cells

    if any(I)
    append!(fy1, FNY.SV[sub2ind(size(FNY), i0[I], j0[I], k0[I])])
    append!(fy2, FNY.SV[sub2ind(size(FNY), i1[I], j0[I], k0[I])])
    append!(fy3, FNY.SV[sub2ind(size(FNY), i0[I], j0[I], k1[I])])
    append!(fy4, FNY.SV[sub2ind(size(FNY), i1[I], j0[I], k1[I])])
    end

    I = (right .== 4) # find  "bigger" cells

    if any(I)
    append!(fy1, FNY.SV[sub2ind(size(FNY), i0[I], j2[I], k0[I])])
    append!(fy2, FNY.SV[sub2ind(size(FNY), i1[I], j2[I], k0[I])])
    append!(fy3, FNY.SV[sub2ind(size(FNY), i0[I], j2[I], k1[I])])
    append!(fy4, FNY.SV[sub2ind(size(FNY), i1[I], j2[I], k1[I])])
    end

    ################################

    I = (front .== 4) # find  "bigger" cells

    if any(I)
    append!(fz1, FNZ.SV[sub2ind(size(FNZ), i0[I], j0[I], k0[I])])
    append!(fz2, FNZ.SV[sub2ind(size(FNZ), i1[I], j0[I], k0[I])])
    append!(fz3, FNZ.SV[sub2ind(size(FNZ), i0[I], j1[I], k0[I])])
    append!(fz4, FNZ.SV[sub2ind(size(FNZ), i1[I], j1[I], k0[I])])
    end

    I = (back .== 4) # find  "bigger" cells

    if any(I)
    append!(fz1, FNZ.SV[sub2ind(size(FNZ), i0[I], j0[I], k2[I])])
    append!(fz2, FNZ.SV[sub2ind(size(FNZ), i1[I], j0[I], k2[I])])
    append!(fz3, FNZ.SV[sub2ind(size(FNZ), i0[I], j1[I], k2[I])])
    append!(fz4, FNZ.SV[sub2ind(size(FNZ), i1[I], j1[I], k2[I])])
    end

    #Convert from n x 1 arrays to vectors
    fx1 = vec(fx1); fx2 = vec(fx2); fx3 = vec(fx3); fx4 = vec(fx4)
    fy1 = vec(fy1); fy2 = vec(fy2); fy3 = vec(fy3); fy4 = vec(fy4)
    fz1 = vec(fz1); fz2 = vec(fz2); fz3 = vec(fz3); fz4 = vec(fz4)


    ###################################################################
    ##
    ## Build interpolation matrix
    ##
    ###################################################################
    #
    # Constraint matrix
    #
    #    f1   f2   f3   f4
    #     I   -I
    #     I        -I
    #     I             -I
    #
    # Null space matrix
    #
    #    f1 I
    #    f2 I
    #    f3 I
    #    f4 I

    nx  = Tn(nnz(FNX)); ny = Tn(nnz(FNY));  nz = Tn(nnz(FNZ));
    Tn1 = one(Tn)
    ## X faces

    n  = Tn(length(fx1)*3)
    i  = vec(repmat([Tn1:n;], 2, 1))
    j  = [vec(repmat(fx1,3,1)); fx2; fx3; fx4]
    v  = [ones(Tf,n); -ones(Tf,n)]
    Cx = sparse(i,j,v,n,nx)

    i  = setdiff([Tn1:nx;], [fx2; fx3; fx4])
    n  = Tn(length(i))
    j  = [Tn1:n;]
    px = zeros(Tn,nx)
    px[i] = j

    k   = indexin(fx1,i)
    j   = [j; vec(repmat(j[k], 3, 1))]
    i   = [i; fx2; fx3; fx4];
    v   = ones(Tf,length(i))
    Nx  = sparse(i,j,v,nx,n)
    
    s   = sum(Nx,1)
    v ./= s[j]
    Qx  = sparse(j,i,v,n,nx)


    ## Y faces
    n  = Tn(length(fy1)*3)
    i  = vec(repmat([Tn1:n;], 2, 1))
    j  = [vec(repmat(fy1,3,1)); fy2; fy3; fy4];
    v  = [ones(Tf,n); -ones(Tf,n)];
    Cy = sparse(i,j,v,n,ny)

    i     = setdiff([Tn1:ny;], [fy2; fy3; fy4]);
    n     = Tn(length(i))
    j     = [Tn1:n;]
    py    = zeros(Tn,ny)
    py[i] = j + maximum(px);

    k   = indexin(fy1,i)
    j   = [j; vec(repmat(j[k], 3, 1))]
    i   = [i; fy2; fy3; fy4]
    v   = ones(Tf,length(i))
    Ny  = sparse(i,j,v,ny,n)
    
    s   = sum(Ny,1)
    v ./= s[j]
    Qy  = sparse(j,i,v,n,ny)

    ## Z faces
    n  = Tn(length(fz1)*3)
    i  = vec(repmat([Tn1:n;], 2, 1))
    j  = [vec(repmat(fz1,3,1)); fz2; fz3; fz4]
    v  = [ones(Tf,n); -ones(Tf,n)]
    Cz = sparse(i,j,v,n,nz)

    i     = setdiff([Tn1:nz;], [fz2; fz3; fz4])
    n     = Tn(length(i))
    j     = [Tn1:n;]
    pz    = zeros(Tn,nz)
    pz[i] = j + maximum(py)

    k   = indexin(fz1,i)
    j   = [j; vec(repmat(j[k], 3, 1))]
    i   = [i; fz2; fz3; fz4]
    v   = ones(Tf,length(i))
    Nz  = sparse(i,j,v,nz,n)
    
    s   = sum(Nz,1)
    v ./= s[j]
    Qz  = sparse(j,i,v,n,nz)

    ## Put it all together

    C  = blkdiag(Cx,Cy,Cz);
    N  = blkdiag(Nx,Ny,Nz);
    Q  = blkdiag(Qx,Qy,Qz);
    p  = [px; py; pz];

    return N,Q,C,p
end

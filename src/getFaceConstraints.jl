export getFaceConstraints

function getFaceConstraints(M::OcTreeMesh)
	if isempty(M.Nf)

		if all(M.S.SV.nzval.==M.S.SV.nzval[1]) # uniform mesh
			M.Nf = speye(sum(M.nf))
			M.Qf = speye(sum(M.nf))
		else
			M.Nf,M.Qf, = getFaceConstraints(M.S)
		end
	end
	return M.Nf,M.Qf
end

function getFaceConstraints(S::SparseArray3D)
  #Eliminate hanging faces
  #   C ... constraint matrix
  #   N ... null space matrix (interpolation: coarse to fine)
  #   Q ... projection matrix (restriction: fine to coarse)
  #                     which satisfy Q * N = I
  #                     (note that Q ~= N')
  #   p ............... lookup table for new face enumeration:
  #                     new face number = p(old face number)
  #                     (returns zero for deleted faces)
  i0,j0,k0,bsz = find3(S);
#  i0       = round(Int64,i0)
#  j0       = round(Int64,j0)
#  k0       = round(Int64,k0)
#  bsz      = round(Int64,bsz)

  i1 = i0 + div(bsz, 2)
  i2 = i0 + bsz
  j1 = j0 + div(bsz, 2)
  j2 = j0 + bsz
  k1 = k0 + div(bsz, 2)
  k2 = k0 + bsz

  upper,lower,left,right,front,back = getNumberOfNeighbors(S);
  nn = [upper; lower; left; right; front; back];
  if ~all( (nn .== 0) | (nn .== 1) | (nn .== 4) )
    error("Implemented only for regularized OcTree meshes")
  end

  FNX,FNY,FNZ = getFaceNumbering(S);

  fx1 = Int64[]; fx2 = Int64[]; fx3 = Int64[]; fx4 = Int64[];
  fy1 = Int64[]; fy2 = Int64[]; fy3 = Int64[]; fy4 = Int64[];
  fz1 = Int64[]; fz2 = Int64[]; fz3 = Int64[]; fz4 = Int64[];


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

  nx = nnz(FNX); ny = nnz(FNY);  nz = nnz(FNZ);

  ## X faces
  n  = length(fx1)*3
  i  = vec(repmat([1:n;], 2, 1))
  j  = [vec(repmat(fx1,3,1)); fx2; fx3; fx4]
  v  = [ones(n); -ones(n)]
  Cx = sparse(i,j,v,n,nx)

  i = setdiff([1:nx;], [fx2; fx3; fx4])
  n = length(i)
  j = [1:n;]
  px = zeros(Int64,nx)
  px[i] = j

  tmp = intersect(fx1,i)
  k   = indexin(tmp,i)
  j   = [j; vec(repmat(j[k], 3, 1))]
  i   = [i; fx2; fx3; fx4];
  v   = ones(length(i))
  Nx  = sparse(i,j,v,nx,n)
  Qx  = spdiagm(vec(1./sum(Nx,1))) * Nx'

  ## Y faces
  n  = length(fy1)*3;
  i  = vec(repmat([1:n;], 2, 1))
  j  = [vec(repmat(fy1,3,1)); fy2; fy3; fy4];
  v  = [ones(n); -ones(n)];
  Cy = sparse(i,j,v,n,ny)

  i     = setdiff([1:ny;], [fy2; fy3; fy4]);
  n     = length(i);
  j     = [1:n;];
  py    = zeros(Int64,ny)
  py[i] = j + maximum(px);

  tmp = intersect(fy1, i)
  k   = indexin(tmp,i)
  j   = [j; vec(repmat(j[k], 3, 1))]
  i   = [i; fy2; fy3; fy4]
  v   = ones(length(i))
  Ny  = sparse(i,j,v,ny,n)
  Qy  = spdiagm(vec(1./sum(Ny,1))) * Ny'

  ## Z faces
  n  = length(fz1)*3
  i  = vec(repmat([1:n;], 2, 1))
  j  = [vec(repmat(fz1,3,1)); fz2; fz3; fz4]
  v  = [ones(n); -ones(n)]
  Cz = sparse(i,j,v,n,nz)

  i     = setdiff([1:nz;], [fz2; fz3; fz4])
  n     = length(i)
  j     = [1:n;]
  pz    = zeros(Int64,nz)
  pz[i] = j + maximum(py)

  tmp = intersect(fz1, i)
  k   = indexin(tmp,i)
  j   = [j; vec(repmat(j[k], 3, 1))]
  i   = [i; fz2; fz3; fz4]
  v   = ones(length(i))
  Nz  = sparse(i,j,v,nz,n)
  Qz  = spdiagm(vec(1./sum(Nz,1))) * Nz'

  ## Put it all together

  C  = blkdiag(Cx,Cy,Cz);
  N  = blkdiag(Nx,Ny,Nz);
  Q  = blkdiag(Qx,Qy,Qz);
  p  = [px; py; pz];

  return N,Q,C,p
end
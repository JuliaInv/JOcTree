export getInterpolationMatrix

function getInterpolationMatrix(M1::OcTreeMesh, M2::OcTreeMesh)

	# check compatibility
	if any(M1.h .!= M2.h)
		error("getInterpolationMatrix requires same cell size")
	end
	if any(M1.x0 .!= M2.x0)
		error("getInterpolationMatrix requires same coordinate origin")
	end
	if any(M1.n .!= M2.n)
		error("getInterpolationMatrix requires same base mesh")
	end

	P = getInterpolationMatrix(M1.S,M2.S)

	return P

end

function getInterpolationMatrix(M1::OcTreeMesh, M2::OcTreeMesh, N::Integer)

	# check compatibility
	if any(M1.h .!= M2.h)
		error("getInterpolationMatrix requires same cell size")
	end
	if any(M1.x0 .!= M2.x0)
		error("getInterpolationMatrix requires same coordinate origin")
	end
	if any(M1.n .!= M2.n)
		error("getInterpolationMatrix requires same base mesh")
	end

	P = getInterpolationMatrix(M1.S,M2.S,N)

	return P

end

getInterpolationMatrix(M1::OcTreeMesh, M2::Array{OcTreeMesh,1}) = map((M) -> getInterpolationMatrix(M1, M), M2)
getInterpolationMatrix(M1::Array{OcTreeMesh,1}, M2::OcTreeMesh) = map((M) -> getInterpolationMatrix(M, M2), M1)

getInterpolationMatrix(M1::OcTreeMesh, M2::Array{OcTreeMesh,1}, N::Integer) = map((M) -> getInterpolationMatrix(M1, M, N), M2)
getInterpolationMatrix(M1::Array{OcTreeMesh,1}, M2::OcTreeMesh, N::Integer) = map((M) -> getInterpolationMatrix(M, M2, N), M1)

function getInterpolationMatrix(S1::SparseArray3D,S2::SparseArray3D)
# P = getInterpolationMatrix(S1::SparseArray3D,S2::SparseArray3D)

C1 = getCellNumbering(S1)
C2 = getCellNumbering(S2)

ii  = Array{Int}(0)
jj  = Array{Int}(0)
val = Array{Float64}(0)

I1 = S1.SV.nzind
I2 = S2.SV.nzind
N1 = falses(prod(S1.sz))
N2 = falses(prod(S2.sz))
N1[I1] = true
N2[I2] = true

# find cells contained both in S1 and S2
i2 = 1
for i1=1:length(I1) # loop over all cells in S1
	while (i2<=length(I2)) && (I1[i1] >= I2[i2]) # loop over all remaining cells in S2
		if (I1[i1] == I2[i2]) && (S1.SV.nzval[i1]==S2.SV.nzval[i2])
			push!(ii, C2.SV.nzval[i2]) # ii = [ii; C2.SV.nzval[i2]]
			push!(jj, C1.SV.nzval[i1])
			push!(val,1.)
		end
		i2 +=1
	end
end


# coarsening of fine cells in S1 to coarser cells in S2
#      S1                                S2
# +-----+--+--+                    +-----------+
# |     |C |F | h/4                |           |
# |  A  +--+--+                    |           |
# |     |D |G | h/4                |           |
# +-----+--+--+       ------->     |     X     | h
# |     |     |                    |           |
# |  B  |  E  | h/2                |           |
# |     |     |                    |           |
# +-----+-----+                    +-----------+
#   h/2   h/2                            h
#
#  Then X = 1/h^2 * [ (h/2)^2*A + (h/2)^2*B + (h/4)^2*C + (h/4)^2*D +
#                     (h/2)^2*E + (h/4)^2*F + (h/4)^2*G  ]

# find cells in S1 not contained in S2
bb = falses(prod(S1.sz))
bb[find(S1.SV .< S2.SV)] = true
bb .|= .!N2
bb .&=  N1
for i1=1:length(I1) # loop over all cells in S1
		if bb[I1[i1]]
			bsz1            = S1.SV.nzval[i1]
			i,j,k           = ind2sub(size(S1),I1[i1])
			bi,bj,bk,bsz2   = findBlocks(S2,i,j,k)  # find parents in S2
	    	push!(ii ,C2.SV[sub2ind(size(C2), bi , bj , bk ),1])
	    	push!(jj ,C1.SV.nzval[i1])
	    	push!(val,(bsz1.^3)./(bsz2.^3))
		end
end

# refinement of coarse cells in S1 to finer cells in S2
#      S2                                S1
# +-----+--+--+                    +-----------+
# |     |C |F | h/4                |           |
# |  A  +--+--+                    |           |
# |     |D |G | h/4                |           |
# +-----+--+--+       <------      |     X     | h
# |     |     |                    |           |
# |  B  |  E  | h/2                |           |
# |     |     |                    |           |
# +-----+-----+                    +-----------+
#   h/2   h/2                            h
#
#  Then INJECT : A=X, B=X, ...., G=X

# find cells in S2 not contained in S1
bb = falses(prod(S1.sz))
bb[find(S1.SV .> S2.SV)] = true
bb .|= .!N1
bb .&=  N2
for i2=1:length(I2); # loop over all cells in S1
		if bb[I2[i2]] # (S1.SV.nzval[i1] < S2.SV.nzval[i2])
			i,j,k      = ind2sub(size(S2),I2[i2])
			bi,bj,bk   = findBlocks(S1,i,j,k) # find parents in S1
			push!(ii ,C2.SV.nzval[i2])
			push!(jj, C1.SV[sub2ind(size(C1) , bi , bj , bk ),1])
			push!(val,1.)
		end
end

return  sparse(vec(ii),vec(jj),vec(val),nnz(S2),nnz(S1))

end

function getInterpolationMatrix(S1::SparseArray3D,S2::SparseArray3D,N::Integer)
# Interpolate cell property on OcTree S1 to cell property on OcTree S2
# using 2^(3*N) point quadrature

n1 = nnz(S1)
n2 = nnz(S2)
# pick faster and less memory consuming version
if 2 * n1 + n2 < n2 * 2 ^ (3 * N)
	return getInterpolationMatrixFilter(S1,S2,N)
else
	return getInterpolationMatrixQuad(S1,S2,N)
end

end


function getInterpolationMatrixFilter(S1::SparseArray3D,S2::SparseArray3D,N::Integer)
# Modification of getInterpolationMatrix(S1::SparseArray3D,S2::SparseArray3D)
# with drop rule for small cells.

n  = 2^N

C1 = getCellNumbering(S1)
C2 = getCellNumbering(S2)

ii  = Array{Int}(0)
jj  = Array{Int}(0)
val = Array{Float64}(0)

I1 = S1.SV.nzind
I2 = S2.SV.nzind
N1 = falses(prod(S1.sz))
N2 = falses(prod(S2.sz))
N1[I1] = true
N2[I2] = true

# find cells contained both in S1 and S2
i2 = 1
for i1=1:length(I1) # loop over all cells in S1
	while (i2<=length(I2)) && (I1[i1] >= I2[i2]) # loop over all remaining cells in S2
		if (I1[i1] == I2[i2]) && (S1.SV.nzval[i1]==S2.SV.nzval[i2])
			push!(ii, C2.SV.nzval[i2]) # ii = [ii; C2.SV.nzval[i2]]
			push!(jj, C1.SV.nzval[i1])
			push!(val,1.)
		end
		i2 +=1
	end
end


# coarsening of fine cells in S1 to coarser cells in S2
#      S1                                S2
# +-----+--+--+                    +-----------+
# |     |C |F | h/4                |           |
# |  A  +--+--+                    |           |
# |     |D |G | h/4                |           |
# +-----+--+--+       ------->     |     X     | h
# |     |     |                    |           |
# |  B  |  E  | h/2                |           |
# |     |     |                    |           |
# +-----+-----+                    +-----------+
#   h/2   h/2                            h
#
#  Then X = 1/h^2 * [ (h/2)^2*A + (h/2)^2*B + (h/4)^2*C + (h/4)^2*D +
#                     (h/2)^2*E + (h/4)^2*F + (h/4)^2*G  ]
#
#  Replace fine cells in S1 if they are smaller than h / 2^N by single cell.

# find cells in S1 not contained in S2
bb = falses(prod(S1.sz))
bb[find(S1.SV .< S2.SV)] = true
bb .|= .!N2
bb .&=  N1
for i1=1:length(I1) # loop over all cells in S1
		if bb[I1[i1]]
			bsz1            = S1.SV.nzval[i1]
			i,j,k           = ind2sub(size(S1),I1[i1])
			bi,bj,bk,bsz2   = findBlocks(S2,i,j,k)  # find parents in S2
			flg,w           = getWeight(i,j,k,bsz1,bi,bj,bk,bsz2,n)
			if flg
		    	push!(ii ,C2.SV[sub2ind(size(C2), bi , bj , bk ),1])
		    	push!(jj ,C1.SV.nzval[i1])
		    	push!(val,w)
			end
		end
end

# refinement of coarse cells in S1 to finer cells in S2
#      S2                                S1
# +-----+--+--+                    +-----------+
# |     |C |F | h/4                |           |
# |  A  +--+--+                    |           |
# |     |D |G | h/4                |           |
# +-----+--+--+       <------      |     X     | h
# |     |     |                    |           |
# |  B  |  E  | h/2                |           |
# |     |     |                    |           |
# +-----+-----+                    +-----------+
#   h/2   h/2                            h
#
#  Then INJECT : A=X, B=X, ...., G=X

# find cells in S2 not contained in S1
bb = falses(prod(S1.sz))
bb[find(S1.SV .> S2.SV)] = true
bb .|= .!N1
bb .&=  N2
for i2=1:length(I2); # loop over all cells in S1
		if bb[I2[i2]] # (S1.SV.nzval[i1] < S2.SV.nzval[i2])
			i,j,k      = ind2sub(size(S2),I2[i2])
			bi,bj,bk   = findBlocks(S1,i,j,k) # find parents in S1
			push!(ii ,C2.SV.nzval[i2])
			push!(jj, C1.SV[sub2ind(size(C1) , bi , bj , bk ),1])
			push!(val,1.)
		end
end

return  sparse(vec(ii),vec(jj),vec(val),nnz(S2),nnz(S1))

end


function getWeight(i1,j1,k1,bsz1,i2,j2,k2,bsz2,n)
# Calculate ratio of volume of cells 1 and 2 where bsz1 < bsz2 and
# filter results based on n^3 quadrature rule on cell 2.
# Note that bsz1, bsz2, n are powers of 2.

if n == 1
	if (i1 == i2) && (j1 == j2) && (k1 == k2)
		return true, 1.0
	else
		return false, 0.0
	end
end

if bsz2 <= n
	return true, (bsz1^3)/(bsz2^3)
end

inc = div(bsz2, n) # integer division is exact since bsz2 and n are powers of 2
if bsz1 < inc
	h2 = div(bsz2, 2)
	#	# sample every inc-th fine cell symmetrically to center of big cell
	#	if (mod(i1-i2,inc) == 0) && (mod(j1-j2,inc) == 0) && (mod(k1-k2,inc) == 0)
	#		return true, (inc^3)/(bsz2^3)
	#	end
	# pick sampling cells symmetric to center of big cell
	if i1 < i2 + h2
		if j1 < j2 + h2
			if k1 < k2 + h2
				if (mod(i1-i2,inc) == 0) && (mod(j1-j2,inc) == 0) && (mod(k1-k2,inc) == 0)
					return true, (inc^3)/(bsz2^3)
				end
			else
				if (mod(i1-i2,inc) == 0) && (mod(j1-j2,inc) == 0) && (mod(k2+bsz2-k1-bsz1,inc) == 0)
					return true, (inc^3)/(bsz2^3)
				end
			end
		else
			if k1 < k2 + h2
				if (mod(i1-i2,inc) == 0) && (mod(j2+bsz2-j1-bsz1,inc) == 0) && (mod(k1-k2,inc) == 0)
					return true, (inc^3)/(bsz2^3)
				end
			else
				if (mod(i1-i2,inc) == 0) && (mod(j2+bsz2-j1-bsz1,inc) == 0) && (mod(k2+bsz2-k1-bsz1,inc) == 0)
					return true, (inc^3)/(bsz2^3)
				end
			end
		end
	else
		if j1 < j2 + h2
			if k1 < k2 + h2
				if (mod(i2+bsz2-i1-bsz1,inc) == 0) && (mod(j1-j2,inc) == 0) && (mod(k1-k2,inc) == 0)
					return true, (inc^3)/(bsz2^3)
				end
			else
				if (mod(i2+bsz2-i1-bsz1,inc) == 0) && (mod(j1-j2,inc) == 0) && (mod(k2+bsz2-k1-bsz1,inc) == 0)
					return true, (inc^3)/(bsz2^3)
				end
			end
		else
			if k1 < k2 + h2
				if (mod(i2+bsz2-i1-bsz1,inc) == 0) && (mod(j2+bsz2-j1-bsz1,inc) == 0) && (mod(k1-k2,inc) == 0)
					return true, (inc^3)/(bsz2^3)
				end
			else
				if (mod(i2+bsz2-i1-bsz1,inc) == 0) && (mod(j2+bsz2-j1-bsz1,inc) == 0) && (mod(k2+bsz2-k1-bsz1,inc) == 0)
					return true, (inc^3)/(bsz2^3)
				end
			end
		end
	end
else
	return true, (bsz1^3)/(bsz2^3)
end

return false, 0.0

end


function getInterpolationMatrixQuad(S1::SparseArray3D,S2::SparseArray3D,N::Integer)
# Straightforward implementation using [2^N]^3 point quadrature in each cell in S2

C1 = getCellNumbering(S1)
C2 = getCellNumbering(S2)

I2,J2,K2,B2 = find3(S2)
n = length(I2)

I = zeros(Int64,0)
J = zeros(Int64,0)
V = zeros(Float64,0)

# iterate over all cells of target OcTree S2
for m = 1:n
	i2 = I2[m]
	j2 = J2[m]
	k2 = K2[m]
	i  = C2.SV[sub2ind(size(C2),i2,j2,k2),1] # row index for P
	d  = getShifts(B2[m], N) # quadrature points relative to (i2,j2,k2)
	v  = 1 / length(d)^3     # quadrature weight (same for all quadrature points)
	# iterate over all quadrature points
	for di = d
		for dj = d
			for dk = d
				i1,j1,k1 = findBlocks(S1, i2+di, j2+dj, k2+dk)
				j = C1.SV[sub2ind(size(C1),i1,j1,k1),1] # column index for P
				push!(I,i)
				push!(J,j)
				push!(V,v)
			end
		end
	end
end

P = sparse(I,J,V,nnz(S2),nnz(S1))

return P

end


function getShifts(bsz::Integer, n::Integer)
# quadrature points relative to lower left front corner of big cell

if n == 0
	# single quadrature point: use lower left front corner
	d = 0
elseif bsz <= 2^n
	# more quadrature points than fine mesh cells: use all fine mesh cells
	d = collect(0:bsz-1)
else
	# less quadrature points than fine mesh cells: sample every inc-th fine cell symmetrically to center of cell
	inc = bsz >> n
	bsz2 = bsz >> 1
	d = [0:inc:bsz2-inc ; bsz-1-(bsz2-inc:-inc:0)]
end

return d

end

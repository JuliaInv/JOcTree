export findBlocks



function findBlocks(S::SparseArray3D,i::Vector,j::Vector,k::Vector)
# [bi,bj,bk,bsz] = findBlocks(S,i,j,k) # search the the blocks in S that covers the indices (i,j,k)
    Tn  = eltype(S.SV.nzval)
	Tn2 = eltype(S.SV.nzind)
	n   = length(i)
	bsz = zeros(Tn,n)
	bi  = copy(i)
	bj  = copy(j)
	bk  = copy(k)

	# find indices outside the domain
	I = find( (bi.<1) .| (S.sz[1].<bi) .| (bj.<1) .| (S.sz[2].<bj) .| (bk.<1) .| (S.sz[3].<bk) )
	if !isempty(I)
		bsz[I] = -one(Tn)
		bi[I]  = -one(Tn2)
		bj[I]  = -one(Tn2)
		bk[I]  = -one(Tn2)
	end

	block = one(Tn2)

	while true
	    Iz = find(bsz.==0)
		if isempty(Iz)
			break
		end
	    t               = S.SV[sub2ind(S.sz,bi[Iz],bj[Iz],bk[Iz])]
	    indnz           = find(t .!= 0)
	    indz            = find(t .== 0)
	    bsz[Iz[indnz]]  = t[indnz]

	    block *= Tn2(2)
		Iz = Iz[indz]
	    bi[Iz] = one(Tn2) .+ block * floor.(Tn2,(i[Iz].-1)/block)
	    bj[Iz] = one(Tn2) .+ block * floor.(Tn2,(j[Iz].-1)/block)
	    bk[Iz] = one(Tn2) .+ block * floor.(Tn2,(k[Iz].-1)/block)
	end
	return bi,bj,bk,bsz
end

function findBlocks(S::SparseArray3D, i::Integer,j::Integer,k::Integer)
	Tn  = eltype(S.SV.nzval)
	Tn2 = eltype(S.SV.nzind)
	if (i < 1) | (i > S.sz[1]) | (j < 1) | (j > S.sz[2]) | (k < 1) | (k > S.sz[3])
		return -one(Tn2), -one(Tn2), -one(Tn2), -one(Tn2)
	end

	bsz = zero(Tn)
	bi  = Tn2(i)
	bj  = Tn2(j)
	bk  = Tn2(k)

	block = one(Tn2)

	while true

		bsz = S.SV[sub2ind(S.sz, bi,bj,bk)]
		if bsz != 0
			break
		end

		block *= Tn2(2)
		bi     = Tn2(1) + block * div(i-1, block)
		bj     = Tn2(1) + block * div(j-1, block)
		bk     = Tn2(1) + block * div(k-1, block)

	end

	return bi,bj,bk,bsz
end

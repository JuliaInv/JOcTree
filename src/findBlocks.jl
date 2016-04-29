export findBlocks



function findBlocks(S::SparseArray3D,i::Vector{Int},j::Vector{Int},k::Vector{Int})
# [bi,bj,bk,bsz] = findBlocks(S,i,j,k) # search the the blocks in S that covers the indices (i,j,k)
	n   = length(i)
	bsz = zeros(Int,n)
	bi  = copy(i)
	bj  = copy(j)
	bk  = copy(k)
	
	# find indices outside the domain
	I = find( (bi.<1) | (S.sz[1].<bi)   |   (bj.<1) | (S.sz[2].<bj) |  (bk.<1) | (S.sz[3].<bk)  )
	if !isempty(I)
		bsz[I] = -1
		bi[I]  = -1
		bj[I]  = -1
		bk[I]  = -1
	end
	
	block = 1
	
	while true
	    Iz = find(bsz.==0)
		if isempty(Iz)
			break 
		end
	    t               = S.SV[sub2ind((S.sz[1],S.sz[2],S.sz[3]),bi[Iz],bj[Iz],bk[Iz])]
	    indnz           = find(t .!= 0)
	    indz            = find(t .== 0)
	    bsz[Iz[indnz]]  = t[indnz]
	    
	    block *= 2
		Iz = Iz[indz]
	    bi[Iz] = 1 .+ block * floor(Integer,(i[Iz].-1)/block)
	    bj[Iz] = 1 .+ block * floor(Integer,(j[Iz].-1)/block)
	    bk[Iz] = 1 .+ block * floor(Integer,(k[Iz].-1)/block)
	end
	return bi,bj,bk,bsz
end

function findBlocks(S,i::Int,j::Int,k::Int)
	
	if (i < 1) | (i > S.sz[1]) | (j < 1) | (j > S.sz[2]) | (k < 1) | (k > S.sz[3]) 
		return -1, -1, -1, -1
	end
	
	bsz = 0
	bi  = copy(i)
	bj  = copy(j)
	bk  = copy(k)
	
	block = 1
	
	while true
		
		bsz = S.SV[sub2ind((S.sz[1],S.sz[2],S.sz[3]),bi,bj,bk),1]
		if bsz != 0
			break
		end
		
	    block *= 2
	    bi     = 1 + block * floor(Integer,(i-1)/block)
	    bj     = 1 + block * floor(Integer,(j-1)/block)
	    bk     = 1 + block * floor(Integer,(k-1)/block)   	      
		
	end
	
	return bi,bj,bk,bsz
end

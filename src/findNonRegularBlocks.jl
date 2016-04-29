export findNonRegularBlocks

function findNonRegularBlocks(S::SparseArray3D)
# findNonRegularBlocks(S)
#

nfac = [
    1/4     1       0;
    1/4     1       1/2;
    3/4     1       0;
    3/4     1       1/2;
    1/4     -1/4    0;
    1/4     -1/4    1/2;
    3/4     -1/4    0;
    3/4     -1/4    1/2;
    1/4     0       -1/4;
    3/4     0       -1/4;
    1/4     1/2     -1/4;
    3/4     1/2     -1/4;
    1/4     0       1;
    3/4     0       1;
    1/4     1/2     1;
    3/4     1/2     1;
    -1/4    0       1/4;
    -1/4    0       3/4;
    -1/4    1/2     1/4;
    -1/4    1/2     3/4;
    1       0       1/4;
    1       0       3/4;
    1       1/2     1/4;
    1       1/2     3/4]


i,j,k,bsz  = find3(S)
meshsize = length(bsz)
isRegular = trues(meshsize)

m1,m2,m3 = S.sz

# check neigbours

for it = 1:meshsize
	bz = bsz[it]
	
	if bz >= 4
		iit = i[it]
		jit = j[it]
		kit = k[it]

		for cnt = 1:24
		   ni = iit + round(Int64, nfac[cnt,1]*bz)
		   nj = jit + round(Int64, nfac[cnt,2]*bz)
		   nk = kit + round(Int64, nfac[cnt,3]*bz)
		
						  

         if 1 <= ni && ni <= m1 &&
            1 <= nj && nj <= m2 &&
            1 <= nk && nk <= m3 
				if S.SV[sub2ind(S.sz,ni,nj,nk),1] != 0
					isRegular[it] = false
					break
			   end
		   end
		   
		end  # cnt

   end # bz >= 4
end  # it

Inr = find(!isRegular)
#Ir  = find(isRegular)


#return i[Inr], j[Inr], k[Inr], bsz[Inr], i[Ir], j[Ir], k[Ir], bsz[Ir]
return Inr, i,j,k, bsz

end # function findNonRegularBlocks

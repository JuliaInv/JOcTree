export getEdgeSize

function getEdgeSize(M::OcTreeMesh)
   if nnz(M.EX) != M.ne[1] || nnz(M.EY) != M.ne[2] || nnz(M.EZ) != M.ne[3]
		M.EX,M.EY,M.EZ = getEdgeSize(M.S)
	end
   return M.EX,M.EY,M.EZ
end


function getEdgeSize(S::SparseArray3D)

m1,m2,m3 = S.sz
I,J,K,BSZ = find3(S)

# allocate
EXi = Int64[]; EXj = Int64[]; EXk = Int64[]; 
EYi = Int64[]; EYj = Int64[]; EYk = Int64[]; 
EZi = Int64[]; EZj = Int64[]; EZk = Int64[]; 
EV  = Int64[];

for bsz = 2.^round(Int64,[0:log2(maximum(BSZ));])
	ind = find(BSZ.==bsz)
	i = I[ind]; j = J[ind]; k = K[ind]
	
	if !isempty(i)
		append!(EV, fill(bsz,4*length(i)) )
		
		# get X edges
		ii  = [  i      ;   i       ;   i       ;    i        ;]
		jj  = [  j      ;   j.+bsz  ;   j       ;   j.+bsz    ;]
		kk  = [  k      ;   k       ;   k.+bsz  ;   k.+bsz    ;]
		append!(EXi, ii)
		append!(EXj, jj)
		append!(EXk, kk)
		
		# get Y edges
		ii    = [  i      ;   i.+bsz   ;   i       ;   i.+bsz   ;]
		jj    = [  j      ;   j        ;   j       ;   j        ;]
		kk    = kk
		append!(EYi, ii)
		append!(EYj, jj)
		append!(EYk, kk)
		
		# get Z edges
		ii    = ii
		jj    = [  j       ;    j       ;   j.+bsz   ;   j.+bsz  ; ]
		kk    = [  k       ;    k       ;   k        ;   k       ; ]
		append!(EZi, ii)
		append!(EZj, jj)
		append!(EZk, kk)
	end
end

EX  = sparse3(EXi,EXj,EXk,EV,[m1;m2+1;m3+1;],(x,y)->x)
EY  = sparse3(EYi,EYj,EYk,EV,[m1+1;m2;m3+1;],(x,y)->x)
EZ  = sparse3(EZi,EZj,EZk,EV,[m1+1;m2+1;m3;],(x,y)->x)


return EX, EY, EZ
end

export createOcTreeFromImage

function createOcTreeFromImage(A::Array{UInt8,3},tol);
# S,A = createOcTreeFromImage(A,tol);
# A - image
# tol - equal intensity color

	m1,m2,m3 = size(A)
	mm = m1*m2*m3

	maxbsz  = minimum(size(A))

# initialize OcTree
	N1   = floor(Integer,m1/maxbsz); N2 = floor(Integer,m2/maxbsz); N3 = floor(Integer,m3/maxbsz)
	ii   = zeros(Int,N1*N2*N3)
	jj   = zeros(Int,N1*N2*N3)
	kk   = zeros(Int,N1*N2*N3)
	bsz  = zeros(Int,N1*N2*N3)
	cnt = 1
	for i=1:maxbsz:m1
		for j=1:maxbsz:m2
			for k = 1:maxbsz:m3
				ii[cnt] = i; jj[cnt] = j; kk[cnt] = k; bsz[cnt] = maxbsz
			 	cnt += 1
			end
		end
	end

	S = Mesh.sparse3(round(Int64,ii),round(Int64,jj),round(Int64,kk),round(Int64,bsz),[m1,m2,m3])
	bszmin = maxbsz
	while true
		println("max blksz = ",maximum(S.SV)," min blksz = ",minimum(nonzeros(S))," number of cells = ",nnz(S))
		nz = nnz(S)
		S = refineOcTreeTol(S,A,tol,bszmin)
		bszmin = floor(Integer,bszmin/2)		
		#S = regularizeOcTree(S)
		if nnz(S) == nz
			return S
		end
	end
	
end


# ------ Refine OcTree
function refineOcTreeTol(S,A,tol,bszmin)

	m1,m2,m3 = S.sz
	ii, jj, kk, bsz = find3(S)
		
	nz = length(ii)
	ti = zeros(Uint,8*nz); tj = zeros(Uint,8*nz); 
	tk = zeros(Uint,8*nz); tb = zeros(Uint,8*nz)
	ti[1:nz] = ii; tj[1:nz] = jj  
	tk[1:nz] = kk; tb[1:nz] = bsz

	cnt = nz+1
	i1 = zero(Uint); j1 = zero(Uint); k1 = zero(Uint); 
	i2 = zero(Uint); j2 = zero(Uint); k2 = zero(Uint)
	for m=1:nz
		if bsz[m] <= bszmin
			i1 = ti[m]; i2 = i1+tb[m]-1
			j1 = tj[m]; j2 = j1+tb[m]-1
			k1 = tk[m]; k2 = k1+tb[m]-1

			#println(i1," ",i2," ",j1," ",j2," ",k1," ",k2," ",tb[m])
			aijk = vec(A[i1:i2,j1:j2,k1:k2])
			blk = copy(tb[m])
			if (maximum(aijk) - minimum(aijk) >= tol)
				tb[m] = blk/2
       
				ti[cnt] = i1+tb[m]
   			tj[cnt] = j1
   			tk[cnt] = k1
				tb[cnt] = blk/2
				cnt += 1
	
				ti[cnt] = i1
   			tj[cnt] = j1+tb[m]
   			tk[cnt] = k1
				tb[cnt] = blk/2
				cnt += 1
	
				ti[cnt] = i1
   			tj[cnt] = j1
   			tk[cnt] = k1+tb[m]
				tb[cnt] = blk/2
				cnt += 1
	
				ti[cnt] = i1+tb[m]
   			tj[cnt] = j1+tb[m]
   			tk[cnt] = k1
				tb[cnt] = blk/2
				cnt += 1

				ti[cnt] = i1+tb[m]
   			tj[cnt] = j1
   			tk[cnt] = k1+tb[m]
				tb[cnt] = blk/2
				cnt += 1
	
				ti[cnt] = i1
   			tj[cnt] = j1+tb[m]
   			tk[cnt] = k1+tb[m]
				tb[cnt] = blk/2
				cnt += 1
	
				ti[cnt] = i1+tb[m]
   			tj[cnt] = j1+tb[m]
   			tk[cnt] = k1+tb[m]
				tb[cnt] = blk/2
				cnt += 1
			end
		end # if
	end # for
	
	ti  = ti[1:cnt-1]
	tj  = tj[1:cnt-1]
	tk  = tk[1:cnt-1]
	tb  = tb[1:cnt-1]
	
	Sr = Mesh.sparse3(round(Int64,ti),round(Int64,tj),round(Int64,tk),round(Int64,tb),[m1,m2,m3])

	ti,tj,tk,tb = find3(Sr)
	
	return Sr

end

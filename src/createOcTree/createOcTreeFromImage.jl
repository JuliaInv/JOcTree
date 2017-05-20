export createOcTreeFromImage

function createOcTreeFromImage(A::Array{UInt8,3},tol);
# S,A = createOcTreeFromImage(A,tol);
# A - image
# tol - equal intensity color

   maxbsz  = minimum(size(A))   
   S = initializeOctree( collect(size(A)) )

	bszmin = maxbsz
	while true
		println("max blksz = ",maximum(S.SV)," min blksz = ",minimum(nonzeros(S))," number of cells = ",nnz(S))
		nz = nnz(S)
		S = refineOcTreeTol(S,A,tol,bszmin)
        bszmin = div(bszmin, 2)
		if nnz(S) == nz
            break
	    end
   end
   return S
   
end


# ------ Refine OcTree
function refineOcTreeTol(S,A,tol,bszmin)

	m1,m2,m3 = S.sz
	ii, jj, kk, bsz = find3(S)
		
	nz = length(ii)
    ti = Array{Int64}(8*nz); tj = Array{Int64}(8*nz)
    tk = Array{Int64}(8*nz); tb = Array{Int64}(8*nz)
	ti[1:nz] = ii; tj[1:nz] = jj  
	tk[1:nz] = kk; tb[1:nz] = bsz
    ii=[] ; jj=[]; kk=[]

	cnt = nz+1
    i1 = zero(Int64); j1 = zero(Int64); k1 = zero(Int64)
    i2 = zero(Int64); j2 = zero(Int64); k2 = zero(Int64)
	for m=1:nz
		if bsz[m] <= bszmin
			i1 = ti[m]; i2 = i1+tb[m]-1
			j1 = tj[m]; j2 = j1+tb[m]-1
			k1 = tk[m]; k2 = k1+tb[m]-1

         minA, maxA = extrema(A[i1:i2,j1:j2,k1:k2])
         blk = tb[m]
         if maxA - minA >= tol
            b2 = div(blk, 2)
            tb[m] = b2
       
            ti[cnt] = i1+b2
   			tj[cnt] = j1
   			tk[cnt] = k1
            tb[cnt] = b2
				cnt += 1
	
				ti[cnt] = i1
            tj[cnt] = j1+b2
   			tk[cnt] = k1
            tb[cnt] = b2
				cnt += 1
	
				ti[cnt] = i1
   			tj[cnt] = j1
            tk[cnt] = k1+b2
            tb[cnt] = b2
				cnt += 1
	
            ti[cnt] = i1+b2
            tj[cnt] = j1+b2
   			tk[cnt] = k1
            tb[cnt] = b2
				cnt += 1

            ti[cnt] = i1+b2
   			tj[cnt] = j1
            tk[cnt] = k1+b2
            tb[cnt] = b2
				cnt += 1
	
				ti[cnt] = i1
            tj[cnt] = j1+b2
            tk[cnt] = k1+b2
            tb[cnt] = b2
				cnt += 1
	
            ti[cnt] = i1+b2
            tj[cnt] = j1+b2
            tk[cnt] = k1+b2
            tb[cnt] = b2
				cnt += 1
			end
		end # if
	end # for
		
   Sr = sparse3(ti[1:cnt-1], tj[1:cnt-1], tk[1:cnt-1], tb[1:cnt-1], [m1,m2,m3])

	return Sr

end

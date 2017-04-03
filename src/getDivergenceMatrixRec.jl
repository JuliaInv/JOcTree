export getDivergenceMatrixRec, getDivergenceMatrix

function getDivergenceMatrix(M::OcTreeMeshFV)
# M.Div = getDivergenceMatrix(M::OcTreeMeshFV) builds face-to-cc divergence operator
        if isempty(M.Div)
			M.Div = getDivergenceMatrixRec(M.S, M.h)
		end
        return M.Div
end

getDivergenceMatrixRec(M) = getDivergenceMatrixRec(M.S,M.h)

function getDivergenceMatrixRec(S,h)
	# [DIV,N,HC,HF,NX,NY,NZ]=getDivergenceMatrixRec(S,h)

	sub2indI(a1, a2, a3, a4) = sub2ind(round(Int64,vec(a1)), round(Int64,vec(a2)), round(Int64,vec(a3)), round(Int64,vec(a4)))

	CN           = getCellNumbering(S)
	FXN,FYN,FZN  = getFaceNumbering(S)
	i,j,k,bsz    = find3(S)
	e1           = ones(length(i))
	upper,lower,left,right,front,back = getNumberOfNeighbors(S)

	#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	#%%% NORMALS ON X-FACES
	#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	ii = sub2indI(CN.sz ,i,j,k) 
	jj = sub2indI(FXN.sz,i,j,k) 
	iind = round(Int64,vec(full(CN.SV[ii,1])))
	jind = round(Int64,vec(full(FXN.SV[jj,1])))
	NX = sparse(iind,jind,-e1,nnz(CN),nnz(FXN))

	I  = find(upper .== 4)
	e  = fill(1,length(I))

	# mark 2nd, 3rd and 4th contributing upper faces
	ii = sub2indI(CN.sz ,i[I],j[I],k[I]) 
	jj = sub2indI(FXN.sz,i[I],j[I]+bsz[I]/2,k[I])
	iind = round(Int64,vec(full(CN.SV[ii,1])))
	jind = round(Int64,vec(full(FXN.SV[jj,1])))
	NX += sparse(iind,jind,-e,nnz(CN),nnz(FXN))

	jj = sub2indI(FXN.sz,i[I],j[I],k[I]+bsz[I]/2)
	jind = round(Int64,vec(full(FXN.SV[jj,1])))
	NX += sparse(iind,jind,-e,nnz(CN),nnz(FXN))

	jj = sub2indI(FXN.sz,i[I],j[I]+bsz[I]/2,k[I]+bsz[I]/2)
	jind = round(Int64,vec(full(FXN.SV[jj,1])))
	NX += sparse(iind,jind,-e,nnz(CN),nnz(FXN))


	# mark 2nd, 3rd and 4th contributing lower faces
	ii = sub2indI(CN.sz ,i,j,k) 
	jj = sub2indI(FXN.sz,i+bsz,j,k) 
	iind = round(Int64,vec(full(CN.SV[ii,1])))
	jind = round(Int64,vec(full(FXN.SV[jj,1])))
	NX += sparse(iind,jind,e1,nnz(CN),nnz(FXN))

	I  = find(lower .== 4)
	e  = fill(1,length(I))

	# mark 2nd, 3rd and 4th contributing upper faces
	ii = sub2indI(CN.sz ,i[I],j[I],k[I]) 
	jj = sub2indI(FXN.sz,i[I]+bsz[I],j[I]+bsz[I]/2,k[I])
	iind = round(Int64,vec(full(CN.SV[ii,1])))
	jind = round(Int64,vec(full(FXN.SV[jj,1])))
	NX += sparse(iind,jind,e,nnz(CN),nnz(FXN))

	jj = sub2indI(FXN.sz,i[I]+bsz[I],j[I],k[I]+bsz[I]/2)
	jind = round(Int64,vec(full(FXN.SV[jj,1])))
	NX += sparse(iind,jind,e,nnz(CN),nnz(FXN))

	jj = sub2indI(FXN.sz,i[I]+bsz[I],j[I]+bsz[I]/2,k[I]+bsz[I]/2)
	jind = round(Int64,vec(full(FXN.SV[jj,1])))
	NX += sparse(iind,jind,e,nnz(CN),nnz(FXN))



	#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	#%%%  NORMALS ON Y-FACES
	#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	ii = sub2indI(CN.sz ,i,j,k) 
	jj = sub2indI(FYN.sz,i,j,k) 
	iind = round(Int64,vec(full(CN.SV[ii,1])))
	jind = round(Int64,vec(full(FYN.SV[jj,1])))
	NY = sparse(iind,jind,-e1,nnz(CN),nnz(FYN))

	jj = sub2indI(FYN.sz,i,j+bsz,k) 
	jind = round(Int64,vec(full(FYN.SV[jj,1])))
	NY += sparse(iind,jind,e1,nnz(CN),nnz(FYN))

	I  = find(left .== 4)
	e  = fill(1,length(I))

	# mark 2nd, 3rd and 4th contributing left faces
	ii = sub2indI(CN.sz ,i[I],j[I],k[I]) 
	jj = sub2indI(FYN.sz,i[I]+bsz[I]/2,j[I],k[I])
	iind = round(Int64,vec(full(CN.SV[ii,1])))
	jind = round(Int64,vec(full(FYN.SV[jj,1])))
	NY += sparse(iind,jind,-e,nnz(CN),nnz(FYN))

	jj = sub2indI(FYN.sz,i[I],j[I],k[I]+bsz[I]/2)
	jind = round(Int64,vec(full(FYN.SV[jj,1])))
	NY += sparse(iind,jind,-e,nnz(CN),nnz(FYN))

	jj = sub2indI(FYN.sz,i[I]+bsz[I]/2,j[I],k[I]+bsz[I]/2)
	jind = round(Int64,vec(full(FYN.SV[jj,1])))
	NY += sparse(iind,jind,-e,nnz(CN),nnz(FYN))


	# mark 2nd, 3rd and 4th contributing right faces
	I  = find(right .== 4)
	e  = fill(1,length(I))

	# mark 2nd, 3rd and 4th contributing upper faces
	ii = sub2indI(CN.sz ,i[I],j[I],k[I]) 
	jj = sub2indI(FYN.sz,i[I]+bsz[I]/2,j[I]+bsz[I],k[I])
	iind = round(Int64,vec(full(CN.SV[ii,1])))
	jind = round(Int64,vec(full(FYN.SV[jj,1])))
	NY += sparse(iind,jind,e,nnz(CN),nnz(FYN))

	jj = sub2indI(FYN.sz,i[I],j[I]+bsz[I],k[I]+bsz[I]/2)
	jind = round(Int64,vec(full(FYN.SV[jj,1])))
	NY += sparse(iind,jind,e,nnz(CN),nnz(FYN))

	jj = sub2indI(FYN.sz,i[I]+bsz[I]/2,j[I]+bsz[I],k[I]+bsz[I]/2)
	jind = round(Int64,vec(full(FYN.SV[jj,1])))
	NY += sparse(iind,jind,e,nnz(CN),nnz(FYN))


	#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	#%%%%  NORMALS ON Z-FACES
	#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	ii = sub2indI(CN.sz ,i,j,k) 
	jj = sub2indI(FZN.sz,i,j,k) 
	iind = round(Int64,vec(full(CN.SV[ii,1])))
	jind = round(Int64,vec(full(FZN.SV[jj,1])))
	NZ = sparse(iind,jind,-e1,nnz(CN),nnz(FZN))

	jj = sub2indI(FZN.sz,i,j,k+bsz) 
	jind = round(Int64,vec(full(FZN.SV[jj,1])))
	NZ += sparse(iind,jind,e1,nnz(CN),nnz(FZN))

	I  = find(front .== 4)
	e  = fill(1,length(I))

	# mark 2nd, 3rd and 4th contributing left faces
	ii = sub2indI(CN.sz ,i[I],j[I],k[I]) 
	jj = sub2indI(FZN.sz,i[I]+bsz[I]/2,j[I],k[I])
	iind = round(Int64,vec(full(CN.SV[ii,1])))
	jind = round(Int64,vec(full(FZN.SV[jj,1])))
	NZ += sparse(iind,jind,-e,nnz(CN),nnz(FZN))

	jj = sub2indI(FZN.sz,i[I],j[I]+bsz[I]/2,k[I])
	jind = round(Int64,vec(full(FZN.SV[jj,1])))
	NZ += sparse(iind,jind,-e,nnz(CN),nnz(FZN))

	jj = sub2indI(FZN.sz,i[I]+bsz[I]/2,j[I]+bsz[I]/2,k[I])
	jind = round(Int64,vec(full(FZN.SV[jj,1])))
	NZ += sparse(iind,jind,-e,nnz(CN),nnz(FZN))


	# mark 2nd, 3rd and 4th contributing right faces
	I  = find(back .== 4)
	e  = fill(1,length(I))

	# mark 2nd, 3rd and 4th contributing upper faces
	ii = sub2indI(CN.sz ,i[I],j[I],k[I]) 
	jj = sub2indI(FZN.sz,i[I]+bsz[I]/2,j[I],k[I]+bsz[I])
	iind = round(Int64,vec(full(CN.SV[ii,1])))
	jind = round(Int64,vec(full(FZN.SV[jj,1])))
	NZ += sparse(iind,jind,e,nnz(CN),nnz(FZN))

	jj = sub2indI(FZN.sz,i[I],j[I]+bsz[I]/2,k[I]+bsz[I])
	jind = round(Int64,vec(full(FZN.SV[jj,1])))
	NZ += sparse(iind,jind,e,nnz(CN),nnz(FZN))

	jj = sub2indI(FZN.sz,i[I]+bsz[I]/2,j[I]+bsz[I]/2,k[I]+bsz[I])
	jind = round(Int64,vec(full(FZN.SV[jj,1])))
	NZ += sparse(iind,jind,e,nnz(CN),nnz(FZN))


	FX,FY,FZ = getFaceSize(S)


	N   = [NX  NY  NZ]

	CSZi = Diagonal(1./(nonzeros(S).^3*prod(h)))

	dd = [nonzeros(FX).^2*h[2]*h[3];
	      nonzeros(FY).^2*h[1]*h[3];
	      nonzeros(FZ).^2*h[1]*h[2];]
	FSZ = Diagonal(dd)

	DIV = CSZi * N * FSZ
	return DIV

end

export getDivergenceMatrixRec, getDivergenceMatrix

"""
    Div = getDivergenceMatrix(M)

    Builds face-to-cc divergence operator

    Input:

        M::OcTreeMeshFV    - The OcTree mesh

    Output:

        Div::SparseMatrixCSC - The discrete divergence matrix
"""
function getDivergenceMatrix(M::OcTreeMeshFV)
# M.Div = getDivergenceMatrix(M::OcTreeMeshFV) builds face-to-cc divergence operator
        if isempty(M.Div)
			M.Div = getDivergenceMatrixRec(M.S, M.h)
		end
        return M.Div
end

getDivergenceMatrixRec(M) = getDivergenceMatrixRec(M.S,M.h)

"""
    Div = getDivergenceMatrixRec(S, h)

    Builds face-to-cc divergence operator

    Input:

        S::SparseArray3D   - Sparse OcTree matrix
        h::Vector{Float64} - Underlying cell size

    Output:

        Div::SparseMatrixCSC - The discrete divergence matrix
"""
function getDivergenceMatrixRec(S::SparseArray3D,h)
#function getDivergenceMatrixRec(M::OcTreeMeshFV)
	# [DIV,N,HC,HF,NX,NY,NZ]=getDivergenceMatrixRec(S,h)

	CN           = getCellNumbering(S)
  FX,FY,FZ, FXN,FYN,FZN = getFaceSizeNumbering(S)
	i,j,k,bsz    = find3(S)
	e1           = ones(length(i))
	upper,lower,left,right,front,back = getNumberOfNeighbors(S)

	#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	#%%% NORMALS ON X-FACES
	#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ii = sub2ind(CN.sz ,i,j,k)

  jj = sub2ind(FXN.sz,i,j,k)
   iind = vec(full(CN.SV[ii]))
   jind = vec(full(FXN.SV[jj]))
	NX = sparse(iind,jind,-e1,nnz(CN),nnz(FXN))

   I  = (upper .== 4)
   e  = ones(Int64,sum(I))

	# mark 2nd, 3rd and 4th contributing upper faces
  jj = sub2ind(FXN.sz,i[I],j[I]+div.(bsz[I],2),k[I])
   iind = vec(full(CN.SV[ii[I]]))
   jind = vec(full(FXN.SV[jj]))
	NX += sparse(iind,jind,-e,nnz(CN),nnz(FXN))

  jj = sub2ind(FXN.sz,i[I],j[I],k[I]+div.(bsz[I],2))
   jind = vec(full(FXN.SV[jj]))
	NX += sparse(iind,jind,-e,nnz(CN),nnz(FXN))

  jj = sub2ind(FXN.sz,i[I],j[I]+div.(bsz[I],2),k[I]+div.(bsz[I],2))
   jind = vec(full(FXN.SV[jj]))
	NX += sparse(iind,jind,-e,nnz(CN),nnz(FXN))


	# mark 2nd, 3rd and 4th contributing lower faces
  jj = sub2ind(FXN.sz,i+bsz,j,k)
   iind = vec(full(CN.SV[ii]))
   jind = vec(full(FXN.SV[jj]))
	NX += sparse(iind,jind,e1,nnz(CN),nnz(FXN))

  I  = (lower .== 4)
  e  = ones(Int64,sum(I))

	# mark 2nd, 3rd and 4th contributing upper faces
  jj = sub2ind(FXN.sz,i[I]+bsz[I],j[I]+div.(bsz[I],2),k[I])
   iind = vec(full(CN.SV[ii[I]]))
   jind = vec(full(FXN.SV[jj]))
	NX += sparse(iind,jind,e,nnz(CN),nnz(FXN))

  jj = sub2ind(FXN.sz,i[I]+bsz[I],j[I],k[I]+div.(bsz[I],2))
   jind = vec(full(FXN.SV[jj]))
	NX += sparse(iind,jind,e,nnz(CN),nnz(FXN))

  jj = sub2ind(FXN.sz,i[I]+bsz[I],j[I]+div.(bsz[I],2),k[I]+div.(bsz[I],2))
   jind = vec(full(FXN.SV[jj]))
	NX += sparse(iind,jind,e,nnz(CN),nnz(FXN))



	#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	#%%%  NORMALS ON Y-FACES
	#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  jj = sub2ind(FYN.sz,i,j,k)
   iind = vec(full(CN.SV[ii]))
   jind = vec(full(FYN.SV[jj]))
	NY = sparse(iind,jind,-e1,nnz(CN),nnz(FYN))

  jj = sub2ind(FYN.sz,i,j+bsz,k) 
   jind = vec(full(FYN.SV[jj]))
	NY += sparse(iind,jind,e1,nnz(CN),nnz(FYN))

  I  = (left .== 4)
  e  = ones(Int64,sum(I))

	# mark 2nd, 3rd and 4th contributing left faces
  jj = sub2ind(FYN.sz,i[I]+div.(bsz[I],2),j[I],k[I])
   iind = vec(full(CN.SV[ii[I]]))
   jind = vec(full(FYN.SV[jj]))
	NY += sparse(iind,jind,-e,nnz(CN),nnz(FYN))

  jj = sub2ind(FYN.sz,i[I],j[I],k[I]+div.(bsz[I],2))
   jind = vec(full(FYN.SV[jj]))
	NY += sparse(iind,jind,-e,nnz(CN),nnz(FYN))

  jj = sub2ind(FYN.sz,i[I]+div.(bsz[I],2),j[I],k[I]+div.(bsz[I],2))
   jind = vec(full(FYN.SV[jj]))
	NY += sparse(iind,jind,-e,nnz(CN),nnz(FYN))


	# mark 2nd, 3rd and 4th contributing right faces
  I  = (right .== 4)
  e  = ones(Int64,sum(I))

	# mark 2nd, 3rd and 4th contributing upper faces
  jj = sub2ind(FYN.sz,i[I]+div.(bsz[I],2),j[I]+bsz[I],k[I])
   iind = vec(full(CN.SV[ii[I]]))
   jind = vec(full(FYN.SV[jj]))
	NY += sparse(iind,jind,e,nnz(CN),nnz(FYN))

  jj = sub2ind(FYN.sz,i[I],j[I]+bsz[I],k[I]+div.(bsz[I],2))
   jind = vec(full(FYN.SV[jj]))
	NY += sparse(iind,jind,e,nnz(CN),nnz(FYN))

  jj = sub2ind(FYN.sz,i[I]+div.(bsz[I],2),j[I]+bsz[I],k[I]+div.(bsz[I],2))
   jind = vec(full(FYN.SV[jj]))
	NY += sparse(iind,jind,e,nnz(CN),nnz(FYN))


	#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	#%%%%  NORMALS ON Z-FACES
	#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  jj = sub2ind(FZN.sz,i,j,k)
   iind = vec(full(CN.SV[ii]))
   jind = vec(full(FZN.SV[jj]))
	NZ = sparse(iind,jind,-e1,nnz(CN),nnz(FZN))

  jj = sub2ind(FZN.sz,i,j,k+bsz)
   jind = vec(full(FZN.SV[jj]))
	NZ += sparse(iind,jind,e1,nnz(CN),nnz(FZN))

  I  = (front .== 4)
  e  = ones(Int64,sum(I))

	# mark 2nd, 3rd and 4th contributing left faces
  jj = sub2ind(FZN.sz,i[I]+div.(bsz[I],2),j[I],k[I])
   iind = vec(full(CN.SV[ii[I]]))
   jind = vec(full(FZN.SV[jj]))
	NZ += sparse(iind,jind,-e,nnz(CN),nnz(FZN))

  jj = sub2ind(FZN.sz,i[I],j[I]+div.(bsz[I],2),k[I])
   jind = vec(full(FZN.SV[jj]))
	NZ += sparse(iind,jind,-e,nnz(CN),nnz(FZN))

  jj = sub2ind(FZN.sz,i[I]+div.(bsz[I],2),j[I]+div.(bsz[I],2),k[I])
   jind = vec(full(FZN.SV[jj]))
	NZ += sparse(iind,jind,-e,nnz(CN),nnz(FZN))


	# mark 2nd, 3rd and 4th contributing right faces
  I  = (back .== 4)
  e  = ones(Int64,sum(I))

	# mark 2nd, 3rd and 4th contributing upper faces
  jj = sub2ind(FZN.sz,i[I]+div.(bsz[I],2),j[I],k[I]+bsz[I])
   iind = vec(full(CN.SV[ii[I]]))
   jind = vec(full(FZN.SV[jj]))
	NZ += sparse(iind,jind,e,nnz(CN),nnz(FZN))

  jj = sub2ind(FZN.sz,i[I],j[I]+div.(bsz[I],2),k[I]+bsz[I])
   jind = vec(full(FZN.SV[jj]))
	NZ += sparse(iind,jind,e,nnz(CN),nnz(FZN))

  jj = sub2ind(FZN.sz,i[I]+div.(bsz[I],2),j[I]+div.(bsz[I],2),k[I]+bsz[I])
   jind = vec(full(FZN.SV[jj]))
	NZ += sparse(iind,jind,e,nnz(CN),nnz(FZN))

  nrow = size(NX,1)
  nc1 = size(NX,2)
  nc2 = size(NY,2)
  nc3 = size(NZ,2)
  N = spzeros(nrow, nc1+nc2+nc3)
  N[:, 1:nc1] = NX
  N[:, nc1+1:nc1+nc2] = NY
  N[:, nc1+nc2+1:nc1+nc2+nc3] = NZ

  CSZi = 1. ./ (nonzeros(S).^3*prod(h))

  FSZ = vcat(nonzeros(FX).^2*(h[2]*h[3]),
             nonzeros(FY).^2*(h[1]*h[3]),
             nonzeros(FZ).^2*(h[1]*h[2]) )

   DIV = DiagTimesMTimesDiag!(CSZi, N, FSZ)

   return DIV
end

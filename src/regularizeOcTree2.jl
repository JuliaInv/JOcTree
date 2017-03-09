export regularizeOcTree2

#using jInv.Mesh
#using Utils.sdiag


function regularizeOcTree2(S::SparseArray3D)
# Refine OcTree to eliminate
# (a) single transitional cells, that is, cells of
#     size bsz which touch cells of both bsz * 2 and bsz / 2;
# (b) cells which are larger than their neighbour by more than a factor
#     of 2 (this covers the function regularizeOcTree).
# Cells are considered neighbours if they share a face, an edge or a
# vertex.


# Regularize. This saves refinement iterations.
S = regularizeOcTree(S)


# Iterate until all cells meet the quality restrictions. Terminates in the
# worst case when all cells are fine cells (bsz = 1).
while true
  
  
  # Unit coefficient node to cell center matrix.
   m1,m2,m3 = S.sz
   i,j,k,bsz = find3(S)
   meshsize = nnz(S)

   ii = vcat( i, i+bsz, i,     i+bsz, i,     i+bsz, i,     i+bsz )
   jj = vcat( j, j,     j+bsz, j+bsz, j,     j,     j+bsz, j+bsz )
   kk = vcat( k, k,     k,     k,     k+bsz, k+bsz, k+bsz, k+bsz )
  
   sizeN = [m1+1, m2+1, m3+1 ]
   ijk = sub2ind( sizeN, ii,jj,kk )
   #[~,~,j] = unique(nind, 'rows')
   j = uniqueidx( ijk )[3]
   ii=0; jj=0; kk=0;
   
   i = vec( repmat((1:meshsize), 8, 1) )
   #N = spones(sparse(i,j,1))
   N = sparse(i, j, ones(Int8,length(i)) )


   # Cell center to cell center matrix via nodal connectivity: cells are
   # neighbours if they share a face, an edge or a vertex.
   #N = spones(N * N')
   # N = spones( A_mul_Bt(N,N) )
   N = A_mul_Bt(N,N)
 #  fill!(N.nzval, 1)
 
  
   # Sparse tabular lookup of cell sizes: the j-th column contains the cell
   # size of the j-th cell and of all its neighbours.
 #  N = sdiag(bsz) * N
  # N = DiagTimesM(bsz, N)
  
   # Convert to logical array: the j-th column is true for all the
   # neighbouring cells who are larger than the j-th cell.
   # Note that in L', the j-th column is true for all the
   # neighbouring cells who are smaller than the j-th cell.
  # L = N .> N'

   #L = speye(Bool, size(N,1), size(N,2))
   m = size(N,1)
   n = size(N,2)
   nnzN = nnz(N)
   Lcolptr = Array{Int64}(length(N.colptr))
   Lrowval = Array{Int64}(nnzN)
   Lnzval  = Array{Bool}(nnzN)
   nnzL = 1
   for ir = 1:meshsize
      j1 = N.colptr[ir]
      j2 = N.colptr[ir+1] - 1
      Lcolptr[ir] = nnzL
      for ic = j1:j2
      	rval = N.rowval[ic]
      	if bsz[rval] > bsz[ir]
      	   Lnzval[ nnzL] = true
      	   Lrowval[nnzL] = rval
      		nnzL += 1
      	end
      end # ic
   end  # ir
   Lcolptr[length(Lcolptr)] = nnzL
   deleteat!(Lnzval,  nnzL:nnzN)
   deleteat!(Lrowval, nnzL:nnzN)
   L = SparseMatrixCSC{Bool,Int64}(m,n,Lcolptr,Lrowval,Lnzval)
   N = 0
   
   
   # Find all cells who have neighbours which are larger and smaller
   #j = full(any(L)) & full(any(L'))
   j = any(L,1) & any(L',1)
  
   # ... and get the neighbours who are larger.
   i,ij = findn(L[:, vec(j)])
  
   # Merge multiple entries.
   i = unique(i)
  
   if isempty(i)
     # Done.
     break
   end
  
   # Refine larger neighbours.
  # S = refineOcTree_(S, tau, 0.9)
   S = splitCells(S, i)
  
end  # while

return S
end # function regularizeOcTree2

#----------------------------------------------------

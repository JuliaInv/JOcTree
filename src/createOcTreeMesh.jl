export initializeOctree, OctreeBox, octreeRegion


function initializeOctree( n::Vector{Int64}) # number of underlying cells
# Initialize octree mesh to the coarsest cells

nb = minimum(n)  # largest cell size
i,j,k = ndgrid(1:nb:n[1], 1:nb:n[2], 1:nb:n[3])

b = fill(nb, length(i))
#S = sparse3(vec(i),vec(j),vec(k), ones(Int,prod(size(i)))*nb, n)
S = sparse3(vec(i),vec(j),vec(k), b, n)

return S	
end # function initializeOctree

#--------------------------------------------

function OctreeBox( S::SparseArray3D,
	                 i1,i2, j1,j2, k1,k2,
	                 cellsize )
# S( i1:i2, j1:j2, k1:k2 ) = cellsize

if  i1>i2 || j1>j2 || k1>k2 ||
	 i1<1 || i2>S.sz[1] ||
	 j1<1 || j2>S.sz[2] ||
	 k1<1 || k2>S.sz[3]
   error("i1>i2 || j1>j2 || k1>k2 ...")
end


while true
   i,j,k,bsz = find3(S)
   tau = zeros(nnz(S))

   for ic = 1:nnz(S)
      ii = i[ic]
      jj = j[ic]
      kk = k[ic]
      bb = bsz[ic]

      if bb <= cellsize
      	continue  # cell is already small
      end


      if ii+bb-1 >= i1 && ii <= i2 &&
         jj+bb-1 >= j1 && jj <= j2 &&
         kk+bb-1 >= k1 && kk <= k2 
         
	      tau[ic] = 1  # large cell to split
      end

   end  # ic
   

   if all(tau .== 0) 
      break  # we are done
   end
    
   S = refineOcTree(S,tau,0.9)

end  # while

return S

end # function OctreeBox

#------------------------------------------------------------

function octreeRegion( S::SparseArray3D,
	                    i1::Vector{Int64}, j1::Vector{Int64}, k1::Vector{Int64},
	                    cellsize )
#  S( i1, j1, k1 ) = cellsize

npts = length(i1)

if length(j1) != npts || length(k1) != npts
   error("length(j1) != npts ...")
end

min_i = minimum(i1)
max_i = maximum(i1)
min_j = minimum(j1)
max_j = maximum(j1)
min_k = minimum(k1)
max_k = maximum(k1)

if  min_i<1 || max_i>S.sz[1] ||
	 min_j<1 || max_j>S.sz[2] ||
	 min_k<1 || max_k>S.sz[3]
   error("min_i<1 || max_i>S.sz[1] ...")
end


while true

   i,j,k,bsz = find3(S)
   tau = zeros(nnz(S))

   for ic = 1:nnz(S)
      ii = i[ic]
      jj = j[ic]
      kk = k[ic]
      bb = bsz[ic]

      if bb <= cellsize
      	continue  # cell is already small
      end

      if ii+bb-1 < min_i || ii > max_i || 
         jj+bb-1 < min_j || jj > max_j || 
         kk+bb-1 < min_k || kk > max_k
      	continue  # cell outside of region
      end


      for ipts = 1:npts

         if ii+bb-1 >= i1[ipts] && ii <= i1[ipts] &&
            jj+bb-1 >= j1[ipts] && jj <= j1[ipts] &&
            kk+bb-1 >= k1[ipts] && kk <= k1[ipts] 

  	         tau[ic] = 1  # large cell to split
  	         break
         end

      end  # ipts

   end  # ic
   

   if all(tau .== 0) 
      break  # we are done
   end
    
   S = refineOcTree(S,tau,0.9)

end  # while


return S

end # function octreeRegion



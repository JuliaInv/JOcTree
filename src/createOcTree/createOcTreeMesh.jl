export initializeOctree, OctreeBox, octreeRegion


function initializeOctree( n::Vector{Int64},   # number of underlying cells
                           bsz::Int64=0)  # cell size (0 for largest)
# Initialize octree mesh to the coarsest cells

if !(ispow2(n[1]) && ispow2(n[2]) && ispow2(n[3]))
   error("n must be power of 2.")
end

if bsz == 0
nb = minimum(n)  # largest cell size
else
   if !ispow2(bsz) || bsz < 0
      error("bsz must be power of 2.")
   elseif any( n .< bsz )
      error("bsz is too large.")
   end
   nb = bsz
end

i,j,k = ndgrid(1:nb:n[1], 1:nb:n[2], 1:nb:n[3])
b = fill(nb, length(i))

S = sparse3(vec(i),vec(j),vec(k), b, n)
return S	
end # function initializeOctree

#--------------------------------------------

function OctreeBox( S::SparseArray3D,
                    i1::Int64, i2::Int64,
                    j1::Int64, j2::Int64,
                    k1::Int64, k2::Int64,
                    cellsize::Int64 )

if  i1>i2 || j1>j2 || k1>k2 ||
	 i1<1 || i2>S.sz[1] ||
	 j1<1 || j2>S.sz[2] ||
	 k1<1 || k2>S.sz[3]
   error("i1>i2 || j1>j2 || k1>k2 ...")
end


while true
   i,j,k,bsz = find3(S)
   nns = nnz(S)
   splitcells = Array{Int64}(nns)
   nsplit = 0  # counter for the number of cells to split.

   for ic = 1:nns
      bb = bsz[ic]

      if bb <= cellsize
      	continue  # cell is already small
      end

      ii = i[ic]
      jj = j[ic]
      kk = k[ic]

      if ii+bb-1 >= i1 && ii <= i2 &&
         jj+bb-1 >= j1 && jj <= j2 &&
         kk+bb-1 >= k1 && kk <= k2 
         
         # large cell to split
         nsplit += 1
         splitcells[nsplit] = ic
      end

   end  # ic
   

   if nsplit == 0 
      break  # we are done
   end
    
   S = splitCells(i,j,k,bsz, S.sz, splitcells[1:nsplit])

end  # while

return S
end # function OctreeBox

#------------------------------------------------------------

function octreeRegion( S::SparseArray3D,
	                    i1::Vector{Int64}, j1::Vector{Int64}, k1::Vector{Int64},
                      cellsize::Int64 )

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

   nns = nnz(S)
   splitcells = Array{Int64}(nns)
   nsplit = 0  # counter for the number of cells to split.

   for ic = 1:nns
      bb = bsz[ic]

      if bb <= cellsize
      	continue  # cell is already small
      end

      ii = i[ic]
      jj = j[ic]
      kk = k[ic]

      if ii+bb-1 < min_i || ii > max_i || 
         jj+bb-1 < min_j || jj > max_j || 
         kk+bb-1 < min_k || kk > max_k
      	continue  # cell outside of region
      end


      for ipts = 1:npts

         if ii+bb-1 >= i1[ipts] && ii <= i1[ipts] &&
            jj+bb-1 >= j1[ipts] && jj <= j1[ipts] &&
            kk+bb-1 >= k1[ipts] && kk <= k1[ipts] 

            # large cell to split
            nsplit += 1
            splitcells[nsplit] = ic
            break
         end

      end  # ipts

   end  # ic
   

   if nsplit == 0 
      break  # we are done
   end
    
   S = splitCells(i,j,k,bsz, S.sz, splitcells[1:nsplit])

end  # while


return S
end # function octreeRegion



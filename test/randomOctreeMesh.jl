
function randomOctreeMesh{Tn<:Integer,Tn2<:Integer}(::Type{Tn}, ::Type{Tn2},
                             n::Tuple, nrand )
# Create random octree mesh
   S = initializeOctree(Tn,Tn2,n)

   ii = rand(1:n[1], nrand)
   jj = rand(1:n[2], nrand)
   kk = rand(1:n[3], nrand)

   S = refineAtPoints(S, ii,jj,kk, 1)
   S = regularizeOcTree(S)

return S
end  # function randomOctreeMesh


function refineAtPoints( S::SparseArray3D,i1::Vector, j1::Vector, k1::Vector,
                         cellsize::Integer )

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

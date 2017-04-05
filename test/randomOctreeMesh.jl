
function randomOctreeMesh( n::Vector, nrand )
# Create random octree mesh
   S = initializeOctree(n)

   ii = rand(1:n[1], nrand)
   jj = rand(1:n[2], nrand)
   kk = rand(1:n[3], nrand)

   S = octreeRegion(S, ii,jj,kk, 1)
   S = regularizeOcTree(S)

return S
end  # function randomOctreeMesh   


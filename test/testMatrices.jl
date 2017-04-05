
using JOcTree
using Base.Test

include("getFunction.jl")
include("randomOctreeMesh.jl") 


# Test the following matrices:
#   getNodalGradientMatrix, getCurlMatrix, getDivergenceMatrix
# These routines (mostly) depend on:
#   getFaceSizeNumbering, getEdgeSizeNumbering, getCellNumbering,
#   getNodalNumbering, getNumberOfNeighbors, sparse3, etc.


println("Testing: getNodalGradientMatrix, getCurlMatrix, getDivergenceMatrix")

#-------------------------------------------

n = [128 , 128 , 128]   # underlying mesh
h = 1 ./ n   # cell size
x0 = [0. , 0. , 0.]

nrand = 5
S = randomOctreeMesh( n, nrand )

M = getOcTreeMeshFV(S, h, x0=x0)

exportUBCOcTreeMesh("mesh.txt", M)

#-------------------------------------------

GRAD = getNodalGradientMatrix(M)
CURL = getCurlMatrix(M)
DIV = getDivergenceMatrix(M)

DIVCURL = DIV * CURL
dropzeros!(DIVCURL)
CURLGRAD = CURL * GRAD
dropzeros!(CURLGRAD)

println()
println("nnz(DIVCURL)   ", nnz(DIVCURL))
println("nnz(CURLGRAD)  ", nnz(CURLGRAD))

@test nnz(DIVCURL) == 0
@test nnz(CURLGRAD) == 0


#-------------------------------------


xyz = getCellCenteredGrid(M)
f = getF( xyz[:,1], xyz[:,2], xyz[:,3] )
exportUBCOcTreeModel("model.txt", M,f)


println()
println("   cells  ||dG||inf  ||dD||inf  ||dC||inf")

ntests = 4
NdifGRAD = Array{Float64}(ntests)
NdifDIV  = Array{Float64}(ntests)
NdifCURL = Array{Float64}(ntests)
ncells = Array{Int64}(ntests)

for i = 1:ntests
   ncells[i] = nnz(S)
   

   GRAD = getNodalGradientMatrix(M)
   xyz = getNodalGrid(M)
   f = getF( xyz[:,1], xyz[:,2], xyz[:,3] )  # on nodes
   Gf = GRAD * f   # Gradient on edges

   # Analytic grad f
   EX,EY,EZ = getEdgeGrids(M)
   dfx = getdFdX( EX[:,1], EX[:,2], EX[:,3] )
   dfy = getdFdY( EY[:,1], EY[:,2], EY[:,3] )
   dfz = getdFdZ( EZ[:,1], EZ[:,2], EZ[:,3] )
   
   aGf = [ dfx; dfy; dfz ]
   difGRAD = Gf - aGf
   NdifGRAD[i] = norm(difGRAD, Inf)

   #-------------------------------------------

   DIV = getDivergenceMatrix(M)
   EX,EY,EZ = getFaceGrids(M)
   fx = getF( EX[:,1], EX[:,2], EX[:,3] )
   fy = getF( EY[:,1], EY[:,2], EY[:,3] )
   fz = getF( EZ[:,1], EZ[:,2], EZ[:,3] )
   
   f = [fx ; fy ; fz]  # on faces
   Df = DIV * f   # Divergence on cell centres

   # Analytic div f
   xyz = getCellCenteredGrid(M)
   dfx = getdFdX( xyz[:,1], xyz[:,2], xyz[:,3] )
   dfy = getdFdY( xyz[:,1], xyz[:,2], xyz[:,3] )
   dfz = getdFdZ( xyz[:,1], xyz[:,2], xyz[:,3] )

   aDf = dfx .+ dfy .+ dfz
   difDIV = Df - aDf
   NdifDIV[i] = norm(difDIV, Inf)

   #-------------------------------------------

   CURL = getCurlMatrix(M)
   EX,EY,EZ = getEdgeGrids(M)
   fx = getF( EX[:,1], EX[:,2], EX[:,3] )
   fy = getF( EY[:,1], EY[:,2], EY[:,3] )
   fz = getF( EZ[:,1], EZ[:,2], EZ[:,3] )
   
   f = [fx ; fy ; fz]  # on edges
   Cf = CURL * f   # Curl on faces

   # Analytic curl f
   EX,EY,EZ = getFaceGrids(M)
   dyz = getdFdZ( EX[:,1], EX[:,2], EX[:,3] )
   dzy = getdFdY( EX[:,1], EX[:,2], EX[:,3] )
   dxz = getdFdZ( EY[:,1], EY[:,2], EY[:,3] )
   dzx = getdFdX( EY[:,1], EY[:,2], EY[:,3] )
   dxy = getdFdY( EZ[:,1], EZ[:,2], EZ[:,3] )
   dyx = getdFdX( EZ[:,1], EZ[:,2], EZ[:,3] )

   aCf = [ -dyz + dzy  ;  dxz - dzx  ;  -dxy + dyx ]
   difCURL = Cf - aCf
   NdifCURL[i] = norm(difCURL, Inf)


   @printf("%8i  %.3e  %.3e  %.3e \n",
               ncells[i], NdifGRAD[i], 
                          NdifDIV[i],  
                          NdifCURL[i] )

   if i < ntests
      S = refineOcTree(S, ones(nnz(S)), 0.1 )
      S = regularizeOcTree(S)
      M = getOcTreeMeshFV(S, h, x0=x0)
   end

end  # i


# Make sure that the inf norm decreases as the number cells increase.

for i = 1 : ntests-1
   @test NdifGRAD[i+1] < NdifGRAD[i]
   @test NdifDIV[i+1]  < NdifDIV[i]
   @test NdifCURL[i+1] < NdifCURL[i]
end  # i


using PyPlot
loglog( ncells, NdifGRAD, "r.-")
loglog( ncells, NdifDIV,  "g.-")
loglog( ncells, NdifCURL, "b.-")

